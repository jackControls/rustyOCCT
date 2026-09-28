//! The Part 21 reader (ISO 10303-21, edition 2): the header, one or more
//! `DATA` sections of simple and complex entity instances, and every
//! parameter kind. It knows no schema: entity and type names are kept as
//! written, strings raw (their `\X\`, `\X2\` and `\X4\` directives are not
//! decoded), and every reference is checked against the instances read.
//! Anything after `END-ISO-10303-21;` (a signature section) is ignored.
use super::StepError;
use std::collections::BTreeMap;

/// Lists may nest this deep; a deeper file is refused rather than risking
/// the stack.
const MAX_DEPTH: usize = 64;

/// One parameter of a record.
#[derive(Debug, Clone, PartialEq)]
pub enum Parameter {
    Integer(i64),
    Real(f64),
    /// The text between the quotes, `''` read as `'`.
    String(String),
    /// `.NAME.`, without the dots.
    Enumeration(String),
    /// `"..."`, the hexadecimal digits.
    Binary(String),
    /// `#n`.
    Reference(u64),
    List(Vec<Parameter>),
    /// `NAME(value)`, such as `LENGTH_MEASURE(1.)`.
    Typed(String, Box<Parameter>),
    /// `$`: an unset optional attribute.
    Unset,
    /// `*`: an attribute redeclared as derived.
    Derived,
}

/// `NAME(parameters)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub name: String,
    pub parameters: Vec<Parameter>,
}

/// An entity instance: one record, or the records of a complex instance in
/// the order written (Part 21 writes them by name).
#[derive(Debug, Clone, PartialEq)]
pub struct Instance {
    /// The line (1-based) where `#n=` stands.
    pub line: usize,
    pub records: Vec<Record>,
    pub complex: bool,
}

impl Instance {
    /// The record of this name, simple or part of a complex instance.
    pub fn record(&self, name: &str) -> Option<&Record> {
        self.records.iter().find(|r| r.name == name)
    }
}

/// A whole exchange structure.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Exchange {
    /// `FILE_DESCRIPTION`, `FILE_NAME`, `FILE_SCHEMA` and any others, in
    /// order.
    pub header: Vec<Record>,
    /// Every instance of every data section by entity number.
    pub instances: BTreeMap<u64, Instance>,
}

impl Exchange {
    /// The schema names of `FILE_SCHEMA`.
    pub fn schemas(&self) -> Vec<&str> {
        self.header
            .iter()
            .filter(|r| r.name == "FILE_SCHEMA")
            .flat_map(|r| r.parameters.first())
            .flat_map(|p| match p {
                Parameter::List(items) => items.iter().collect(),
                other => vec![other],
            })
            .filter_map(|p| match p {
                Parameter::String(s) => Some(s.as_str()),
                _ => None,
            })
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Keyword(String),
    Name(u64),
    Integer(i64),
    Real(f64),
    Str(String),
    Enumeration(String),
    Binary(String),
    Open,
    Close,
    Comma,
    Semicolon,
    Equals,
    Dollar,
    Star,
}

struct Lexer<'a> {
    text: &'a [u8],
    at: usize,
    line: usize,
}

fn syntax(line: usize, what: &'static str) -> StepError {
    StepError::Syntax { line, what }
}

impl Lexer<'_> {
    fn error(&self, what: &'static str) -> StepError {
        syntax(self.line, what)
    }

    /// Whitespace and comments.
    fn skip(&mut self) -> Result<(), StepError> {
        loop {
            match self.text.get(self.at) {
                Some(b'\n') => {
                    self.line += 1;
                    self.at += 1;
                }
                Some(c) if *c <= b' ' => self.at += 1,
                Some(b'/') if self.text.get(self.at + 1) == Some(&b'*') => {
                    let start = self.line;
                    self.at += 2;
                    loop {
                        match self.text.get(self.at) {
                            None => return Err(syntax(start, "the end of a comment")),
                            Some(b'*') if self.text.get(self.at + 1) == Some(&b'/') => {
                                self.at += 2;
                                break;
                            }
                            Some(c) => {
                                if *c == b'\n' {
                                    self.line += 1;
                                }
                                self.at += 1;
                            }
                        }
                    }
                }
                _ => return Ok(()),
            }
        }
    }

    fn take_while(&mut self, f: impl Fn(u8) -> bool) -> &str {
        let start = self.at;
        while self.text.get(self.at).is_some_and(|c| f(*c)) {
            self.at += 1;
        }
        // Only ASCII bytes pass every predicate used here.
        std::str::from_utf8(&self.text[start..self.at]).unwrap_or("")
    }

    /// The next token and the line it starts on, or `None` at the end.
    fn next(&mut self) -> Result<Option<(usize, Token)>, StepError> {
        self.skip()?;
        let line = self.line;
        let Some(&c) = self.text.get(self.at) else {
            return Ok(None);
        };
        let token = match c {
            b'(' | b')' | b',' | b';' | b'=' | b'$' | b'*' => {
                self.at += 1;
                match c {
                    b'(' => Token::Open,
                    b')' => Token::Close,
                    b',' => Token::Comma,
                    b';' => Token::Semicolon,
                    b'=' => Token::Equals,
                    b'$' => Token::Dollar,
                    _ => Token::Star,
                }
            }
            b'#' => {
                self.at += 1;
                let digits = self.take_while(|c| c.is_ascii_digit());
                if digits.is_empty() {
                    return Err(self.error("an entity number after #"));
                }
                Token::Name(
                    digits
                        .parse()
                        .map_err(|_| syntax(line, "an entity number within 64 bits"))?,
                )
            }
            b'\'' => {
                self.at += 1;
                let mut bytes = Vec::new();
                loop {
                    match self.text.get(self.at) {
                        None => return Err(syntax(line, "the end of a string")),
                        Some(b'\'') if self.text.get(self.at + 1) == Some(&b'\'') => {
                            bytes.push(b'\'');
                            self.at += 2;
                        }
                        Some(b'\'') => {
                            self.at += 1;
                            break;
                        }
                        // `\\` is a backslash, and `\S\` takes the next
                        // character whatever it is, an apostrophe included;
                        // both are kept raw.
                        Some(b'\\') => {
                            let n = match (self.text.get(self.at + 1), self.text.get(self.at + 2)) {
                                (Some(b'\\'), _) => 2,
                                (Some(b'S'), Some(b'\\')) if self.at + 3 < self.text.len() => 4,
                                _ => 1,
                            };
                            let raw = &self.text[self.at..self.at + n];
                            self.line += raw.iter().filter(|c| **c == b'\n').count();
                            bytes.extend_from_slice(raw);
                            self.at += n;
                        }
                        Some(&c) => {
                            if c == b'\n' {
                                self.line += 1;
                            }
                            bytes.push(c);
                            self.at += 1;
                        }
                    }
                }
                Token::Str(String::from_utf8_lossy(&bytes).into_owned())
            }
            b'"' => {
                self.at += 1;
                let digits = self.take_while(|c| c.is_ascii_hexdigit()).to_string();
                if self.text.get(self.at) != Some(&b'"') {
                    return Err(self.error("the end of a binary"));
                }
                self.at += 1;
                Token::Binary(digits)
            }
            b'.' => {
                self.at += 1;
                let name = self
                    .take_while(|c| c.is_ascii_alphanumeric() || c == b'_')
                    .to_ascii_uppercase();
                if name.is_empty() || self.text.get(self.at) != Some(&b'.') {
                    return Err(self.error("an enumeration"));
                }
                self.at += 1;
                Token::Enumeration(name)
            }
            b'+' | b'-' | b'0'..=b'9' => self.number(line)?,
            c if c.is_ascii_alphabetic() || c == b'_' || c == b'!' => {
                let start = self.at;
                self.at += 1;
                self.take_while(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-');
                let word = std::str::from_utf8(&self.text[start..self.at])
                    .unwrap_or("")
                    .to_ascii_uppercase();
                // Hyphens belong only to the file's delimiters.
                if word.contains('-') && word != "ISO-10303-21" && word != "END-ISO-10303-21" {
                    return Err(syntax(line, "a keyword"));
                }
                Token::Keyword(word)
            }
            _ => return Err(self.error("a token")),
        };
        Ok(Some((line, token)))
    }

    /// `[+-]digits[.digits][E[+-]digits]`: a real when it has a point or
    /// an exponent, else an integer.
    fn number(&mut self, line: usize) -> Result<Token, StepError> {
        let start = self.at;
        if matches!(self.text.get(self.at), Some(b'+' | b'-')) {
            self.at += 1;
        }
        if self.take_while(|c| c.is_ascii_digit()).is_empty() {
            return Err(syntax(line, "digits in a number"));
        }
        let mut real = false;
        if self.text.get(self.at) == Some(&b'.') {
            real = true;
            self.at += 1;
            self.take_while(|c| c.is_ascii_digit());
        }
        if matches!(self.text.get(self.at), Some(b'E' | b'e')) {
            real = true;
            self.at += 1;
            if matches!(self.text.get(self.at), Some(b'+' | b'-')) {
                self.at += 1;
            }
            if self.take_while(|c| c.is_ascii_digit()).is_empty() {
                return Err(syntax(line, "digits in an exponent"));
            }
        }
        let text = std::str::from_utf8(&self.text[start..self.at]).unwrap_or("");
        if real {
            match text.parse::<f64>() {
                Ok(x) if x.is_finite() => Ok(Token::Real(x)),
                _ => Err(syntax(line, "a finite real")),
            }
        } else {
            text.parse()
                .map(Token::Integer)
                .map_err(|_| syntax(line, "an integer within 64 bits"))
        }
    }
}

struct Parser<'a> {
    lexer: Lexer<'a>,
    peeked: Option<(usize, Token)>,
}

impl Parser<'_> {
    fn peek(&mut self) -> Result<Option<&Token>, StepError> {
        if self.peeked.is_none() {
            self.peeked = self.lexer.next()?;
        }
        Ok(self.peeked.as_ref().map(|(_, t)| t))
    }

    fn take(&mut self, what: &'static str) -> Result<(usize, Token), StepError> {
        self.peek()?;
        self.peeked
            .take()
            .ok_or_else(|| syntax(self.lexer.line, what))
    }

    fn expect(&mut self, token: Token, what: &'static str) -> Result<usize, StepError> {
        let (line, found) = self.take(what)?;
        if found == token {
            Ok(line)
        } else {
            Err(syntax(line, what))
        }
    }

    fn keyword(&mut self, word: &str, what: &'static str) -> Result<(), StepError> {
        self.expect(Token::Keyword(word.into()), what).map(|_| ())
    }

    fn parameter(&mut self, depth: usize) -> Result<Parameter, StepError> {
        let (line, token) = self.take("a parameter")?;
        Ok(match token {
            Token::Integer(i) => Parameter::Integer(i),
            Token::Real(x) => Parameter::Real(x),
            Token::Str(s) => Parameter::String(s),
            Token::Enumeration(e) => Parameter::Enumeration(e),
            Token::Binary(b) => Parameter::Binary(b),
            Token::Name(n) => Parameter::Reference(n),
            Token::Dollar => Parameter::Unset,
            Token::Star => Parameter::Derived,
            Token::Open => Parameter::List(self.list(depth + 1, line)?),
            Token::Keyword(name) => {
                self.expect(Token::Open, "( after a type name")?;
                if depth + 1 > MAX_DEPTH {
                    return Err(syntax(line, "a nesting depth within 64"));
                }
                let value = self.parameter(depth + 1)?;
                self.expect(Token::Close, ") after a typed parameter")?;
                Parameter::Typed(name, Box::new(value))
            }
            _ => return Err(syntax(line, "a parameter")),
        })
    }

    /// The items of a list whose `(` was just read, and its `)`.
    fn list(&mut self, depth: usize, line: usize) -> Result<Vec<Parameter>, StepError> {
        if depth > MAX_DEPTH {
            return Err(syntax(line, "a nesting depth within 64"));
        }
        let mut items = Vec::new();
        if self.peek()? == Some(&Token::Close) {
            self.take(")")?;
            return Ok(items);
        }
        loop {
            items.push(self.parameter(depth)?);
            match self.take(", or )")? {
                (_, Token::Comma) => continue,
                (_, Token::Close) => return Ok(items),
                (line, _) => return Err(syntax(line, ", or )")),
            }
        }
    }

    /// `NAME(parameters)`.
    fn record(&mut self) -> Result<Record, StepError> {
        let (line, token) = self.take("an entity name")?;
        let Token::Keyword(name) = token else {
            return Err(syntax(line, "an entity name"));
        };
        let open = self.expect(Token::Open, "( after an entity name")?;
        let parameters = self.list(1, open)?;
        Ok(Record { name, parameters })
    }
}

/// Every reference in a parameter.
fn references(p: &Parameter, out: &mut Vec<u64>) {
    match p {
        Parameter::Reference(n) => out.push(*n),
        Parameter::List(items) => items.iter().for_each(|i| references(i, out)),
        Parameter::Typed(_, value) => references(value, out),
        _ => {}
    }
}

/// Read a Part 21 exchange structure.
pub fn read(text: &[u8]) -> Result<Exchange, StepError> {
    let mut p = Parser {
        lexer: Lexer {
            text,
            at: 0,
            line: 1,
        },
        peeked: None,
    };
    p.keyword("ISO-10303-21", "ISO-10303-21")?;
    p.expect(Token::Semicolon, "; after ISO-10303-21")?;
    p.keyword("HEADER", "HEADER")?;
    p.expect(Token::Semicolon, "; after HEADER")?;
    let mut exchange = Exchange::default();
    while p.peek()? != Some(&Token::Keyword("ENDSEC".into())) {
        exchange.header.push(p.record()?);
        p.expect(Token::Semicolon, "; after a header entity")?;
    }
    p.keyword("ENDSEC", "ENDSEC")?;
    p.expect(Token::Semicolon, "; after ENDSEC")?;
    loop {
        match p.take("DATA or END-ISO-10303-21")? {
            (_, Token::Keyword(k)) if k == "END-ISO-10303-21" => {
                p.expect(Token::Semicolon, "; after END-ISO-10303-21")?;
                break;
            }
            (_, Token::Keyword(k)) if k == "DATA" => {}
            (line, _) => return Err(syntax(line, "DATA or END-ISO-10303-21")),
        }
        // Edition 3 names a section and its schemas: read and ignored.
        if p.peek()? == Some(&Token::Open) {
            let (line, _) = p.take("(")?;
            p.list(1, line)?;
        }
        p.expect(Token::Semicolon, "; after DATA")?;
        loop {
            let (line, token) = p.take("an entity instance or ENDSEC")?;
            let number = match token {
                Token::Name(n) => n,
                Token::Keyword(k) if k == "ENDSEC" => break,
                _ => return Err(syntax(line, "an entity instance or ENDSEC")),
            };
            p.expect(Token::Equals, "= after an entity number")?;
            let complex = p.peek()? == Some(&Token::Open);
            let records = if complex {
                p.take("(")?;
                let mut records = Vec::new();
                while p.peek()? != Some(&Token::Close) {
                    records.push(p.record()?);
                }
                p.take(")")?;
                if records.is_empty() {
                    return Err(syntax(line, "a record in a complex instance"));
                }
                records
            } else {
                vec![p.record()?]
            };
            p.expect(Token::Semicolon, "; after an entity instance")?;
            let instance = Instance {
                line,
                records,
                complex,
            };
            if exchange.instances.insert(number, instance).is_some() {
                return Err(syntax(line, "a unique entity number"));
            }
        }
        p.expect(Token::Semicolon, "; after ENDSEC")?;
    }
    let mut targets = Vec::new();
    for (number, instance) in &exchange.instances {
        targets.clear();
        for r in &instance.records {
            r.parameters
                .iter()
                .for_each(|q| references(q, &mut targets));
        }
        if let Some(missing) = targets.iter().find(|t| !exchange.instances.contains_key(t)) {
            return Err(StepError::Reference {
                entity: *number,
                target: *missing,
            });
        }
    }
    Ok(exchange)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn file(data: &str) -> String {
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('','',(''),(''),'','','');\n\
             FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{data}ENDSEC;\nEND-ISO-10303-21;\n"
        )
    }

    #[test]
    fn every_parameter_kind() {
        let text = file(
            "#1=A(1,-2,+3,1.,-0.5E-3,2.5e2,'it''s','',.T.,.UNSPECIFIED.,\"0FF\",#2,$,*,(),(1,(2,#2)),\
             LENGTH_MEASURE(1.E-07),T((#1)));\n#2=B();\n",
        );
        let x = read(text.as_bytes()).unwrap();
        assert_eq!(x.schemas(), vec!["AUTOMOTIVE_DESIGN"]);
        let r = &x.instances[&1].records[0];
        use Parameter::*;
        assert_eq!(
            r.parameters,
            vec![
                Integer(1),
                Integer(-2),
                Integer(3),
                Real(1.0),
                Real(-0.5e-3),
                Real(250.0),
                String("it's".into()),
                String(std::string::String::new()),
                Enumeration("T".into()),
                Enumeration("UNSPECIFIED".into()),
                Binary("0FF".into()),
                Reference(2),
                Unset,
                Derived,
                List(vec![]),
                List(vec![Integer(1), List(vec![Integer(2), Reference(2)])]),
                Typed("LENGTH_MEASURE".into(), Box::new(Real(1e-7))),
                Typed("T".into(), Box::new(List(vec![Reference(1)]))),
            ]
        );
        assert!(x.instances[&2].records[0].parameters.is_empty());
    }

    #[test]
    fn string_directives_stay_raw() {
        let text = file("#1=A('A\\S\\'h','x\\\\','\\X2\\00E9\\X0\\');\n");
        let x = read(text.as_bytes()).unwrap();
        assert_eq!(
            x.instances[&1].records[0].parameters,
            [
                Parameter::String("A\\S\\'h".into()),
                Parameter::String("x\\\\".into()),
                Parameter::String("\\X2\\00E9\\X0\\".into()),
            ]
        );
    }

    #[test]
    fn complex_instances_comments_and_order() {
        let text = file(
            "/* a comment\n over two lines */\n#20=( A(1) B() C(#10) );\n\n#10 = D ( 'x\ny' ) ;\n",
        );
        let x = read(text.as_bytes()).unwrap();
        let c = &x.instances[&20];
        assert!(c.complex);
        assert_eq!(c.line, 10);
        assert_eq!(
            c.records
                .iter()
                .map(|r| r.name.as_str())
                .collect::<Vec<_>>(),
            ["A", "B", "C"]
        );
        assert_eq!(
            c.record("C").unwrap().parameters,
            [Parameter::Reference(10)]
        );
        assert_eq!(x.instances[&10].line, 12);
        assert_eq!(
            x.instances[&10].records[0].parameters,
            [Parameter::String("x\ny".into())]
        );
    }

    #[test]
    fn several_sections_and_a_signature() {
        let text = "ISO-10303-21;HEADER;FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));ENDSEC;\
                    DATA;#1=A(#2);ENDSEC;DATA(('S',('X')));#2=B();ENDSEC;END-ISO-10303-21;\
                    SIGNATURE ignored";
        let x = read(text.as_bytes()).unwrap();
        assert_eq!(x.instances.len(), 2);
        assert_eq!(x.schemas(), vec!["CONFIG_CONTROL_DESIGN"]);
    }

    #[test]
    fn errors_name_the_line_and_expectation() {
        let bad = |data: &str| read(file(data).as_bytes()).unwrap_err();
        assert_eq!(
            bad("#1=A(#2);\n"),
            StepError::Reference {
                entity: 1,
                target: 2
            }
        );
        assert_eq!(
            bad("#1=A();\n#1=B();\n"),
            syntax(9, "a unique entity number")
        );
        assert_eq!(bad("#1=A(1E999);\n"), syntax(8, "a finite real"));
        assert_eq!(
            bad("#1=A(99999999999999999999);\n"),
            syntax(8, "an integer within 64 bits")
        );
        assert_eq!(bad("#1=A('open);\n"), syntax(8, "the end of a string"));
        assert_eq!(bad("#1=A(.T);\n"), syntax(8, "an enumeration"));
        assert_eq!(bad("#1=A(1 2);\n"), syntax(8, ", or )"));
        assert_eq!(bad("#1=();\n"), syntax(8, "a record in a complex instance"));
        assert_eq!(bad("/* open\n"), syntax(8, "the end of a comment"));
        let deep = format!("#1=A({}{});\n", "(".repeat(70), ")".repeat(70));
        assert_eq!(bad(&deep), syntax(8, "a nesting depth within 64"));
        assert!(matches!(read(b""), Err(StepError::Syntax { line: 1, .. })));
        assert!(matches!(
            read(b"ISO-10303-21;HEADER;ENDSEC;DATA;"),
            Err(StepError::Syntax { .. })
        ));
    }

    #[test]
    fn every_fixture_reads() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/step");
        let mut names: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        names.sort();
        assert_eq!(names.len(), 29);
        for path in names {
            let x = read(&std::fs::read(&path).unwrap()).unwrap();
            assert!(x.instances.len() > 20, "{}", path.display());
            assert_eq!(x.schemas().len(), 1);
        }
    }
}
