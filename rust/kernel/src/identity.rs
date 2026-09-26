//! Value identity for topological entities (see `IDENTITY_AND_HISTORY.md`).
//!
//! An [`EntityId`] is the 128-bit FNV-1a digest of a [`Derivation`]: the
//! operation that created the entity, its role, a canonical ordinal and its
//! parents (caller labels, profile indices or other entity ids). Ids are
//! plain values: no pointer, counter or binary64 coordinate enters them, so
//! the same operation on the same labelled inputs yields the same ids on any
//! platform and in any process.
//!
//! Version-1 encoding (integers little-endian): `b"RSID"`, version `1`,
//! operation `u64`, operation kind `u8`, entity kind `u8`, role `u8`, ordinal
//! `u32`, parent count `u32`, then each parent: `1, label u64`;
//! `2, boundary u32, element u8, index u32`; or `3, 16 id bytes`. The id is
//! the FNV-1a-128 digest of these bytes, stored big-endian. The encoding and
//! digest are pinned by `fixtures/identity-vectors.tsv`, produced by the
//! independent `tools/identity_reference.py`.
use crate::{Error, Result};
use std::fmt;
use std::str::FromStr;

const MAGIC: &[u8; 4] = b"RSID";
const VERSION: u8 = 1;
const FNV_OFFSET: u128 = 0x6c62272e07bb014262b821756295c58d;
const FNV_PRIME: u128 = 0x0000000001000000000000000000013B;

/// 128-bit digest of a [`Derivation`]. Displays as 32 lowercase hex digits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId(pub [u8; 16]);

/// A caller label on a profile boundary, segment or vertex. Labels are the
/// roots of derivations: edits that keep labels keep ids.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InputLabel(pub u64);

/// The version of an operation's behaviour (H8). A history records the
/// level it ran at; changing an operation's behaviour requires a new level,
/// and the previous level stays callable. Levels never enter a derivation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AlgorithmLevel(pub u32);

impl AlgorithmLevel {
    /// The first recorded level, introduced with the cell-complex model; no
    /// earlier behaviour is replayable.
    pub const FIRST: Self = Self(1);
    /// The level new operations run at.
    pub const CURRENT: Self = Self::FIRST;
    /// Every level this build replays, oldest first.
    pub const REPLAYABLE: &'static [Self] = &[Self::FIRST];
}

/// A caller-supplied operation identity, e.g. an application feature id.
/// Convenience constructors without one use [`OperationId::UNSPECIFIED`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct OperationId(pub u64);

impl OperationId {
    pub const UNSPECIFIED: Self = Self(0);
}

/// Encoded codes are append-only; never renumber a released variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperationKind {
    Extrude,
    Transform,
    /// Entities supplied through `Topology::from_parts`.
    External,
    /// A history composed with `History::then`; never part of a derivation.
    Composite,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntityKind {
    Vertex,
    Edge,
    Face,
    Body,
    /// A bounded region of space (TOPOLOGY_MODEL.md).
    Region,
}

impl EntityKind {
    /// Topological dimension: vertex 0, edge 1, face 2, body 3.
    pub fn dimension(self) -> u8 {
        match self {
            Self::Vertex => 0,
            Self::Edge => 1,
            Self::Face => 2,
            Self::Body | Self::Region => 3,
        }
    }
}

/// Why an operation created an entity. "Bottom" and "top" mean the start and
/// end sides of an extrusion, whichever offset is lower.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Role {
    StartCap,
    EndCap,
    Wall,
    BottomEdge,
    TopEdge,
    Vertical,
    /// Retired by the seamless cell model (TOPOLOGY_MODEL.md, T1); its code
    /// stays reserved so encodings never change meaning.
    Seam,
    BottomVertex,
    TopVertex,
    /// Retired with `Seam`; its code stays reserved.
    SeamVertex,
    Body,
    External,
    /// A bounded region, generated from every boundary of its profile.
    Region,
}

/// A profile boundary, one of its segments, or one of its vertices, by
/// stored (counter-clockwise) index. Used as a parent when a boundary has no
/// labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProfileElement {
    Boundary,
    Segment(u32),
    Vertex(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Parent {
    Label(InputLabel),
    Profile {
        boundary: u32,
        element: ProfileElement,
    },
    Entity(EntityId),
}

impl Parent {
    /// Dimension of the parent, for history dimension checks: a boundary
    /// stands for the planar region it encloses, a segment is a curve, a
    /// vertex a point. Entity and label parents are resolved elsewhere.
    pub fn profile_dimension(self) -> Option<u8> {
        match self {
            Self::Profile {
                element: ProfileElement::Vertex(_),
                ..
            } => Some(0),
            Self::Profile {
                element: ProfileElement::Segment(_),
                ..
            } => Some(1),
            Self::Profile {
                element: ProfileElement::Boundary,
                ..
            } => Some(2),
            _ => None,
        }
    }
}

/// Why an entity has its id. Stored beside the id, never recomputed from
/// geometry.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Derivation {
    pub operation: OperationId,
    pub kind: OperationKind,
    pub entity: EntityKind,
    pub role: Role,
    /// Exact canonical ordinal among siblings with the same role and parents.
    pub ordinal: u32,
    /// Canonical order, defined per operation.
    pub parents: Vec<Parent>,
}

fn operation_code(kind: OperationKind) -> u8 {
    match kind {
        OperationKind::Extrude => 1,
        OperationKind::Transform => 2,
        OperationKind::External => 3,
        OperationKind::Composite => 4,
    }
}

fn entity_code(kind: EntityKind) -> u8 {
    match kind {
        EntityKind::Vertex => 1,
        EntityKind::Edge => 2,
        EntityKind::Face => 3,
        EntityKind::Body => 4,
        EntityKind::Region => 5,
    }
}

pub(crate) fn role_code(role: Role) -> u8 {
    match role {
        Role::StartCap => 1,
        Role::EndCap => 2,
        Role::Wall => 3,
        Role::BottomEdge => 4,
        Role::TopEdge => 5,
        Role::Vertical => 6,
        Role::Seam => 7,
        Role::BottomVertex => 8,
        Role::TopVertex => 9,
        Role::SeamVertex => 10,
        Role::Body => 11,
        Role::External => 12,
        Role::Region => 13,
    }
}

const OPERATIONS: [OperationKind; 4] = [
    OperationKind::Extrude,
    OperationKind::Transform,
    OperationKind::External,
    OperationKind::Composite,
];
const ENTITIES: [EntityKind; 5] = [
    EntityKind::Vertex,
    EntityKind::Edge,
    EntityKind::Face,
    EntityKind::Body,
    EntityKind::Region,
];
const ROLES: [Role; 13] = [
    Role::StartCap,
    Role::EndCap,
    Role::Wall,
    Role::BottomEdge,
    Role::TopEdge,
    Role::Vertical,
    Role::Seam,
    Role::BottomVertex,
    Role::TopVertex,
    Role::SeamVertex,
    Role::Body,
    Role::External,
    Role::Region,
];

/// Append one parent's encoding.
pub(crate) fn encode_parent(out: &mut Vec<u8>, parent: &Parent) {
    match parent {
        Parent::Label(label) => {
            out.push(1);
            out.extend_from_slice(&label.0.to_le_bytes());
        }
        Parent::Profile { boundary, element } => {
            out.push(2);
            out.extend_from_slice(&boundary.to_le_bytes());
            let (tag, index) = match element {
                ProfileElement::Boundary => (0u8, 0u32),
                ProfileElement::Segment(i) => (1, *i),
                ProfileElement::Vertex(i) => (2, *i),
            };
            out.push(tag);
            out.extend_from_slice(&index.to_le_bytes());
        }
        Parent::Entity(id) => {
            out.push(3);
            out.extend_from_slice(&id.0);
        }
    }
}

/// Parent count, then each parent.
pub(crate) fn encode_parents(out: &mut Vec<u8>, parents: &[Parent]) {
    out.extend_from_slice(&(parents.len() as u32).to_le_bytes());
    for parent in parents {
        encode_parent(out, parent);
    }
}

impl Derivation {
    /// The version-1 byte encoding described in the module documentation.
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(32 + 16 * self.parents.len());
        out.extend_from_slice(MAGIC);
        out.push(VERSION);
        out.extend_from_slice(&self.operation.0.to_le_bytes());
        out.push(operation_code(self.kind));
        out.push(entity_code(self.entity));
        out.push(role_code(self.role));
        out.extend_from_slice(&self.ordinal.to_le_bytes());
        encode_parents(&mut out, &self.parents);
        out
    }

    /// Inverse of [`Derivation::encode`]; rejects any other byte string.
    pub fn decode(bytes: &[u8]) -> Result<Self> {
        let bad = Error::InvalidTopology("malformed derivation encoding");
        let mut r = Reader(bytes);
        if r.take(4).ok_or(bad.clone())? != MAGIC || r.u8().ok_or(bad.clone())? != VERSION {
            return Err(bad);
        }
        let operation = OperationId(r.u64().ok_or(bad.clone())?);
        let pick = |code: u8, n: usize| (code as usize).checked_sub(1).filter(|i| *i < n);
        let kind =
            OPERATIONS[pick(r.u8().ok_or(bad.clone())?, OPERATIONS.len()).ok_or(bad.clone())?];
        let entity =
            ENTITIES[pick(r.u8().ok_or(bad.clone())?, ENTITIES.len()).ok_or(bad.clone())?];
        let role = ROLES[pick(r.u8().ok_or(bad.clone())?, ROLES.len()).ok_or(bad.clone())?];
        let ordinal = r.u32().ok_or(bad.clone())?;
        let count = r.u32().ok_or(bad.clone())?;
        let mut parents = Vec::new();
        for _ in 0..count {
            let parent = match r.u8().ok_or(bad.clone())? {
                1 => Parent::Label(InputLabel(r.u64().ok_or(bad.clone())?)),
                2 => {
                    let boundary = r.u32().ok_or(bad.clone())?;
                    let tag = r.u8().ok_or(bad.clone())?;
                    let index = r.u32().ok_or(bad.clone())?;
                    let element = match (tag, index) {
                        (0, 0) => ProfileElement::Boundary,
                        (1, i) => ProfileElement::Segment(i),
                        (2, i) => ProfileElement::Vertex(i),
                        _ => return Err(bad),
                    };
                    Parent::Profile { boundary, element }
                }
                3 => {
                    let mut id = [0u8; 16];
                    id.copy_from_slice(r.take(16).ok_or(bad.clone())?);
                    Parent::Entity(EntityId(id))
                }
                _ => return Err(bad),
            };
            parents.push(parent);
        }
        if !r.0.is_empty() {
            return Err(bad);
        }
        Ok(Self {
            operation,
            kind,
            entity,
            role,
            ordinal,
            parents,
        })
    }

    pub fn id(&self) -> EntityId {
        EntityId(fnv1a128(&self.encode()))
    }
}

struct Reader<'a>(&'a [u8]);
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        if self.0.len() < n {
            return None;
        }
        let (head, tail) = self.0.split_at(n);
        self.0 = tail;
        Some(head)
    }
    fn u8(&mut self) -> Option<u8> {
        self.take(1).map(|b| b[0])
    }
    fn u32(&mut self) -> Option<u32> {
        self.take(4)
            .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
    }
    fn u64(&mut self) -> Option<u64> {
        self.take(8)
            .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
    }
}

/// FNV-1a with a 128-bit state, big-endian output. Stable across Rust
/// versions and platforms, unlike the standard library's hashers.
pub fn fnv1a128(bytes: &[u8]) -> [u8; 16] {
    let mut hash = FNV_OFFSET;
    for byte in bytes {
        hash ^= u128::from(*byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash.to_be_bytes()
}

impl fmt::Display for EntityId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for byte in self.0 {
            write!(f, "{byte:02x}")?;
        }
        Ok(())
    }
}

impl FromStr for EntityId {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        let bad = Error::InvalidTopology("malformed entity id");
        if s.len() != 32
            || !s
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
        {
            return Err(bad);
        }
        let mut id = [0u8; 16];
        for (i, byte) in id.iter_mut().enumerate() {
            *byte = u8::from_str_radix(&s[2 * i..2 * i + 2], 16).map_err(|_| bad.clone())?;
        }
        Ok(Self(id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unhex(s: &str) -> Vec<u8> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect()
    }

    #[test]
    fn fixed_vectors_pin_the_encoding_and_digest() {
        let rows = include_str!("../../fixtures/identity-vectors.tsv");
        let mut count = 0;
        for line in rows.lines().filter(|l| !l.starts_with('#')) {
            let fields: Vec<&str> = line.split('\t').collect();
            let (name, encoding, id) = (fields[0], unhex(fields[1]), fields[2]);
            let derivation = Derivation::decode(&encoding).unwrap_or_else(|_| panic!("{name}"));
            assert_eq!(derivation.encode(), encoding, "{name}: round trip");
            assert_eq!(derivation.id().to_string(), id, "{name}: digest");
            assert_eq!(id.parse::<EntityId>().unwrap(), derivation.id(), "{name}");
            count += 1;
        }
        assert_eq!(count, 12);
        // Published FNV-1a-128 test vectors.
        assert_eq!(
            EntityId(fnv1a128(b"")).to_string(),
            "6c62272e07bb014262b821756295c58d"
        );
        assert_eq!(
            EntityId(fnv1a128(b"a")).to_string(),
            "d228cb696f1a8caf78912b704e4a8964"
        );
        assert_eq!(
            EntityId(fnv1a128(b"foobar")).to_string(),
            "343e1662793c64bf6f0d3597ba446f18"
        );
    }

    #[test]
    fn malformed_encodings_and_ids_are_rejected() {
        let good = Derivation {
            operation: OperationId(3),
            kind: OperationKind::Extrude,
            entity: EntityKind::Edge,
            role: Role::Vertical,
            ordinal: 0,
            parents: vec![Parent::Profile {
                boundary: 1,
                element: ProfileElement::Vertex(2),
            }],
        }
        .encode();
        assert!(Derivation::decode(&good).is_ok());
        for cut in 0..good.len() {
            assert!(Derivation::decode(&good[..cut]).is_err(), "prefix {cut}");
        }
        let mut extra = good.clone();
        extra.push(0);
        assert!(Derivation::decode(&extra).is_err());
        for (index, value) in [
            (0, b'X'),
            (4, 2),
            (13, 0),
            (13, 5),
            (14, 6),
            (15, 14),
            (24, 9),
            (29, 3),
        ] {
            let mut bad = good.clone();
            bad[index] = value;
            assert!(Derivation::decode(&bad).is_err(), "byte {index} = {value}");
        }
        for s in ["", "00", &"g".repeat(32), &"A".repeat(32), &"0".repeat(33)] {
            assert!(s.parse::<EntityId>().is_err(), "{s}");
        }
    }
}
