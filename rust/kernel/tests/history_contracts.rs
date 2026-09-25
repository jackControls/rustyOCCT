//! Contract checkers and composition laws on hand-written histories (M0).
use rusty_occt::attributes::{
    self, Attribute, AttributeKey, AttributeMap, AttributeOutcome, OnMerge, OnModify, OnSplit,
    OnTransform, OutcomeResult, Policy, PolicyTable,
};
use rusty_occt::history::{
    check, EntityInfo, EntitySet, Geometry, History, HistoryIssueKind as K, Relation, Resolution,
};
use rusty_occt::identity::{
    EntityId, EntityKind, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use rusty_occt::topology::{Curve3, Orientation, Surface};
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;

fn id(n: u8) -> EntityId {
    EntityId([n; 16])
}
fn p(x: f64) -> Point3 {
    Point3::new(x, 0.0, 0.0)
}
fn line(a: f64, b: f64) -> Geometry {
    Geometry::Curve(Curve3::LineSegment {
        start: p(a),
        end: p(b),
    })
}
fn plane(z: f64) -> Geometry {
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, z),
        Vec3::Z,
        Vec3::X,
        Tolerance::default(),
    )
    .unwrap();
    Geometry::Surface {
        surface: Surface::Plane(frame),
        orientation: Orientation::Forward,
    }
}
fn cylinder() -> Geometry {
    Geometry::Surface {
        surface: Surface::Cylinder {
            frame: Frame3::xy(),
            radius: 1.0,
        },
        orientation: Orientation::Forward,
    }
}
fn info(kind: EntityKind, ordinal: u32, geometry: Geometry) -> EntityInfo {
    EntityInfo {
        kind,
        ordinal,
        geometry,
        structure: Vec::new(),
    }
}
fn set(body: u8, entities: Vec<(u8, EntityInfo)>) -> EntitySet {
    EntitySet {
        body: id(body),
        tolerance: Tolerance::default(),
        entities: entities.into_iter().map(|(n, e)| (id(n), e)).collect(),
    }
}
fn history(kind: OperationKind, from: u8, to: u8, relations: Vec<Relation>) -> History {
    History::new(
        OperationId(7),
        kind,
        vec![id(from)],
        vec![id(to)],
        relations,
        vec![],
    )
}
fn kinds(issues: &[rusty_occt::history::HistoryIssue]) -> Vec<(K, EntityId)> {
    issues.iter().map(|i| (i.kind, i.id)).collect()
}

/// Input body 100: edge 1 on x in [0, 2], face 2 at z=0, vertex 3, face 4.
fn input() -> EntitySet {
    set(
        100,
        vec![
            (1, info(EntityKind::Edge, 0, line(0.0, 2.0))),
            (2, info(EntityKind::Face, 0, plane(0.0))),
            (3, info(EntityKind::Vertex, 0, Geometry::Point(p(5.0)))),
            (4, info(EntityKind::Face, 0, plane(1.0))),
        ],
    )
}

/// Output body 101: edge 1 split into 10 and 11, face 2 unchanged, vertex 3
/// deleted, face 4 modified, and edge 12 generated from face 2.
fn output() -> EntitySet {
    set(
        101,
        vec![
            (10, info(EntityKind::Edge, 0, line(0.0, 1.0))),
            (11, info(EntityKind::Edge, 1, line(1.0, 2.0))),
            (2, info(EntityKind::Face, 0, plane(0.0))),
            (4, info(EntityKind::Face, 0, plane(1.0))),
            (12, info(EntityKind::Edge, 0, line(7.0, 8.0))),
        ],
    )
}

fn clean() -> Vec<Relation> {
    vec![
        Relation::Split {
            from: id(1),
            into: vec![id(10), id(11)],
        },
        Relation::Unchanged { id: id(2) },
        Relation::Deleted { id: id(3) },
        Relation::Modified {
            from: id(4),
            to: id(4),
        },
        Relation::Generated {
            from: vec![Parent::Entity(id(2))],
            to: id(12),
            role: Role::Vertical,
        },
    ]
}

#[test]
fn a_clean_history_has_no_issues_and_resolves_every_input() {
    let h = history(OperationKind::Extrude, 100, 101, clean());
    assert_eq!(check(&[input()], &[output()], &h), vec![]);
    assert_eq!(h.resolve(id(1)), Resolution::Split(vec![id(10), id(11)]));
    assert_eq!(h.resolve(id(2)), Resolution::Same(id(2)));
    assert_eq!(h.resolve(id(3)), Resolution::Deleted);
    assert_eq!(h.resolve(id(4)), Resolution::Same(id(4)));
    assert_eq!(h.resolve(id(99)), Resolution::Unknown);
    // Relations are canonically ordered by construction.
    let mut shuffled = clean();
    shuffled.reverse();
    assert_eq!(history(OperationKind::Extrude, 100, 101, shuffled), h);
}

#[test]
fn every_issue_kind_is_reported_on_its_offending_id() {
    let run = |relations: Vec<Relation>, input: EntitySet, output: EntitySet| {
        kinds(&check(
            &[input],
            &[output],
            &history(OperationKind::Extrude, 100, 101, relations),
        ))
    };
    let without = |n: usize| {
        let mut r = clean();
        r.remove(n);
        r
    };
    let with = |extra: Relation| {
        let mut r = clean();
        r.push(extra);
        r
    };
    let replace = |n: usize, new: Relation| {
        let mut r = clean();
        r[n] = new;
        r
    };
    // Dropping the split loses a source and both targets.
    assert_eq!(
        run(without(0), input(), output()),
        vec![
            (K::MissingSource, id(1)),
            (K::MissingTarget, id(10)),
            (K::MissingTarget, id(11))
        ]
    );
    assert_eq!(
        run(with(Relation::Deleted { id: id(2) }), input(), output()),
        vec![(K::DuplicateSource, id(2)), (K::DeletedStillPresent, id(2))]
    );
    assert_eq!(
        run(with(Relation::Unchanged { id: id(50) }), input(), output()),
        vec![(K::DanglingId, id(50))]
    );
    // Split pieces must stay on the parent's line, in ordinal order.
    let mut moved = output();
    moved.entities.get_mut(&id(11)).unwrap().geometry = line(1.0, 3.0);
    assert_eq!(run(clean(), input(), moved.clone()), vec![]);
    moved.entities.get_mut(&id(11)).unwrap().geometry = Geometry::Curve(Curve3::LineSegment {
        start: p(1.0),
        end: Point3::new(2.0, 1e-3, 0.0),
    });
    assert_eq!(
        run(clean(), input(), moved),
        vec![(K::SplitSupportDiffers, id(11))]
    );
    assert_eq!(
        run(
            replace(
                0,
                Relation::Split {
                    from: id(1),
                    into: vec![id(11), id(10)]
                }
            ),
            input(),
            output()
        ),
        vec![
            (K::OrdinalNotCanonical, id(10)),
            (K::OrdinalNotCanonical, id(11))
        ]
    );
    let mut single = output();
    single.entities.remove(&id(11));
    assert_eq!(
        run(
            replace(
                0,
                Relation::Split {
                    from: id(1),
                    into: vec![id(10)]
                }
            ),
            input(),
            single
        ),
        vec![(K::InvalidArity, id(1))]
    );
    let mut changed = output();
    changed.entities.get_mut(&id(2)).unwrap().geometry = plane(0.5);
    assert_eq!(
        run(clean(), input(), changed),
        vec![(K::UnchangedGeometryDiffers, id(2))]
    );
    let mut restructured = output();
    restructured.entities.get_mut(&id(2)).unwrap().structure =
        vec![vec![(id(10), Orientation::Forward)]];
    assert_eq!(
        run(clean(), input(), restructured),
        vec![(K::UnchangedGeometryDiffers, id(2))]
    );
    let mut family = output();
    family.entities.get_mut(&id(4)).unwrap().geometry = cylinder();
    assert_eq!(
        run(clean(), input(), family),
        vec![(K::ModifiedFamilyDiffers, id(4))]
    );
    let mut renamed = output();
    let face = renamed.entities.remove(&id(4)).unwrap();
    renamed.entities.insert(id(40), face);
    assert_eq!(
        run(
            replace(
                3,
                Relation::Modified {
                    from: id(4),
                    to: id(40)
                }
            ),
            input(),
            renamed
        ),
        vec![(K::ModifiedIdChanged, id(4))]
    );
    // A face cannot be split into an edge; a vertex cannot generate a face.
    let mut faces = output();
    faces.entities.get_mut(&id(11)).unwrap().kind = EntityKind::Face;
    faces.entities.get_mut(&id(11)).unwrap().geometry = plane(3.0);
    assert_eq!(
        run(clean(), input(), faces),
        vec![(K::DimensionMismatch, id(11))]
    );
    let mut swept = output();
    swept.entities.get_mut(&id(12)).unwrap().kind = EntityKind::Face;
    assert_eq!(
        run(
            replace(
                4,
                Relation::Generated {
                    from: vec![Parent::Profile {
                        boundary: 0,
                        element: ProfileElement::Vertex(0)
                    }],
                    to: id(12),
                    role: Role::Wall,
                }
            ),
            input(),
            swept
        ),
        vec![(K::DimensionMismatch, id(12))]
    );
    assert_eq!(
        run(
            replace(
                4,
                Relation::Generated {
                    from: vec![],
                    to: id(12),
                    role: Role::Wall
                }
            ),
            input(),
            output()
        ),
        vec![(K::InvalidArity, id(12))]
    );
    // Merged pieces must lie on the result's support.
    let merged_out = set(101, vec![(20, info(EntityKind::Edge, 0, line(0.0, 3.0)))]);
    let merged_in = set(
        100,
        vec![
            (21, info(EntityKind::Edge, 0, line(0.0, 1.0))),
            (
                22,
                info(
                    EntityKind::Edge,
                    1,
                    Geometry::Curve(Curve3::LineSegment {
                        start: Point3::new(1.0, 1.0, 0.0),
                        end: p(3.0),
                    }),
                ),
            ),
        ],
    );
    let merge = vec![Relation::Merged {
        from: vec![id(21), id(22)],
        into: id(20),
    }];
    assert_eq!(
        run(merge, merged_in, merged_out),
        vec![(K::MergedSupportDiffers, id(22))]
    );
    // Bodies, duplicate ids and relation order.
    assert_eq!(
        kinds(&check(
            &[input()],
            &[output()],
            &history(OperationKind::Extrude, 100, 102, clean())
        )),
        vec![(K::BodyMismatch, id(101)), (K::BodyMismatch, id(102))]
    );
    let mut twin = output();
    twin.body = id(103);
    twin.entities.retain(|k, _| *k == id(2));
    let mut two = history(OperationKind::Extrude, 100, 101, clean());
    two.output_bodies.push(id(103));
    assert_eq!(
        kinds(&check(&[input()], &[output(), twin], &two)),
        vec![(K::DuplicateId, id(2))]
    );
    let mut unsorted = history(OperationKind::Extrude, 100, 101, clean());
    unsorted.relations.reverse();
    assert!(kinds(&check(&[input()], &[output()], &unsorted))
        .iter()
        .all(|(k, _)| *k == K::RelationOrder));
    assert!(!check(&[input()], &[output()], &unsorted).is_empty());
}

// ------------------------------------------------------------------ composition

/// Edge chain: body a has edge 1 on [0, 4]; step one splits it at 2, step two
/// splits the first half at 1 and deletes nothing, step three merges all.
fn bodies() -> [EntitySet; 4] {
    let edge =
        |n: u8, a: f64, b: f64, ordinal: u32| (n, info(EntityKind::Edge, ordinal, line(a, b)));
    [
        set(
            1,
            vec![
                edge(1, 0.0, 4.0, 0),
                (9, info(EntityKind::Face, 0, plane(0.0))),
            ],
        ),
        set(
            2,
            vec![
                edge(2, 0.0, 2.0, 0),
                edge(3, 2.0, 4.0, 1),
                (9, info(EntityKind::Face, 0, plane(0.0))),
            ],
        ),
        set(
            3,
            vec![
                edge(4, 0.0, 1.0, 0),
                edge(5, 1.0, 2.0, 1),
                edge(3, 2.0, 4.0, 1),
                (9, info(EntityKind::Face, 0, plane(0.0))),
                (6, info(EntityKind::Vertex, 0, Geometry::Point(p(1.0)))),
            ],
        ),
        set(
            4,
            vec![
                edge(7, 0.0, 4.0, 0),
                (9, info(EntityKind::Face, 0, plane(0.0))),
            ],
        ),
    ]
}

fn steps() -> [History; 3] {
    [
        History::new(
            OperationId(1),
            OperationKind::Extrude,
            vec![id(1)],
            vec![id(2)],
            vec![
                Relation::Split {
                    from: id(1),
                    into: vec![id(2), id(3)],
                },
                Relation::Unchanged { id: id(9) },
            ],
            vec![],
        ),
        History::new(
            OperationId(2),
            OperationKind::Extrude,
            vec![id(2)],
            vec![id(3)],
            vec![
                Relation::Split {
                    from: id(2),
                    into: vec![id(4), id(5)],
                },
                Relation::Unchanged { id: id(3) },
                Relation::Unchanged { id: id(9) },
                Relation::Generated {
                    from: vec![Parent::Entity(id(2))],
                    to: id(6),
                    role: Role::Vertical,
                },
            ],
            vec![],
        ),
        History::new(
            OperationId(3),
            OperationKind::Extrude,
            vec![id(3)],
            vec![id(4)],
            vec![
                Relation::Merged {
                    from: vec![id(3), id(4), id(5)],
                    into: id(7),
                },
                Relation::Unchanged { id: id(9) },
                Relation::Deleted { id: id(6) },
            ],
            vec![],
        ),
    ]
}

#[test]
fn each_step_and_every_composition_checks_clean() {
    let b = bodies();
    let s = steps();
    for (k, h) in s.iter().enumerate() {
        assert_eq!(check(&b[k..=k], &b[k + 1..=k + 1], h), vec![], "step {k}");
    }
    let ab = s[0].then(&s[1]).unwrap();
    let bc = s[1].then(&s[2]).unwrap();
    let abc = ab.then(&s[2]).unwrap();
    assert_eq!(abc, s[0].then(&bc).unwrap(), "associativity");
    assert_eq!(check(&b[0..1], &b[2..3], &ab), vec![]);
    assert_eq!(check(&b[1..2], &b[3..4], &bc), vec![]);
    assert_eq!(check(&b[0..1], &b[3..4], &abc), vec![]);
    // The split of edge 1 flattens through the second split.
    assert_eq!(
        ab.resolve(id(1)),
        Resolution::Split(vec![id(4), id(5), id(3)])
    );
    // The generated vertex is reported against the original input.
    assert!(ab.relations.contains(&Relation::Generated {
        from: vec![Parent::Entity(id(1))],
        to: id(6),
        role: Role::Vertical
    }));
    // Split then merge of the same pieces is Modified, never Unchanged.
    assert_eq!(abc.resolve(id(1)), Resolution::Same(id(7)));
    assert!(abc.relations.contains(&Relation::Modified {
        from: id(1),
        to: id(7)
    }));
    assert_eq!(abc.resolve(id(9)), Resolution::Same(id(9)));
    assert!(abc.relations.contains(&Relation::Unchanged { id: id(9) }));
    assert_eq!(abc.kind, OperationKind::Composite);
    // Merge with other inputs reports every parent, and deletion composes.
    assert_eq!(
        bc.resolve(id(2)),
        Resolution::Merged {
            into: id(7),
            with: vec![id(3)]
        }
    );
    assert_eq!(
        bc.resolve(id(3)),
        Resolution::Merged {
            into: id(7),
            with: vec![id(2)]
        }
    );
}

#[test]
fn identity_steps_are_neutral_and_mismatched_chains_are_rejected() {
    let s = steps();
    let keep = |body: u8, ids: &[u8]| {
        History::new(
            OperationId(9),
            OperationKind::Transform,
            vec![id(body)],
            vec![id(body)],
            ids.iter()
                .map(|n| Relation::Modified {
                    from: id(*n),
                    to: id(*n),
                })
                .collect(),
            vec![],
        )
    };
    let before = keep(1, &[1, 9]).then(&s[0]).unwrap();
    let after = s[0].then(&keep(2, &[2, 3, 9])).unwrap();
    // A transform's Modified turns Unchanged into Modified; relations agree otherwise.
    for h in [&before, &after] {
        assert_eq!(h.resolve(id(1)), Resolution::Split(vec![id(2), id(3)]));
        assert_eq!(h.resolve(id(9)), Resolution::Same(id(9)));
        assert!(h.relations.contains(&Relation::Modified {
            from: id(9),
            to: id(9)
        }));
    }
    let invalid = |r: Result<History, rusty_occt::history::HistoryIssue>| r.unwrap_err().kind;
    assert_eq!(invalid(s[0].then(&s[2])), K::CompositionInvalid);
    let mut incomplete = s[1].clone();
    incomplete
        .relations
        .retain(|r| *r != Relation::Unchanged { id: id(3) });
    assert_eq!(invalid(s[0].then(&incomplete)), K::CompositionInvalid);
    // A split whose piece merges with another input is many-to-many.
    let mixed = History::new(
        OperationId(4),
        OperationKind::Extrude,
        vec![id(2)],
        vec![id(3)],
        vec![
            Relation::Merged {
                from: vec![id(3), id(9)],
                into: id(8),
            },
            Relation::Unchanged { id: id(2) },
        ],
        vec![],
    );
    assert_eq!(invalid(s[0].then(&mixed)), K::CompositionInvalid);
}

// ------------------------------------------------------------------ attributes

fn attr(key: u64, value: &[u8]) -> Attribute {
    Attribute {
        key: AttributeKey(key),
        value: value.to_vec(),
    }
}

fn policy(modify: OnModify, transform: OnTransform, split: OnSplit, merge: OnMerge) -> Policy {
    Policy {
        on_modify: modify,
        on_transform: transform,
        on_split: split,
        on_merge: merge,
    }
}

#[test]
fn attribute_outcomes_follow_every_policy_and_relation_kind() {
    use OutcomeResult::*;
    // Input: 1 split, 2 unchanged, 3 deleted, 4 modified, 5 and 6 merged.
    let relations = vec![
        Relation::Split {
            from: id(1),
            into: vec![id(10), id(11)],
        },
        Relation::Unchanged { id: id(2) },
        Relation::Deleted { id: id(3) },
        Relation::Modified {
            from: id(4),
            to: id(4),
        },
        Relation::Merged {
            from: vec![id(5), id(6)],
            into: id(12),
        },
    ];
    for kind in [OperationKind::Extrude, OperationKind::Transform] {
        for modify in [OnModify::Keep, OnModify::Drop] {
            for transform in [OnTransform::Keep, OnTransform::Drop, OnTransform::Recompute] {
                for split in [OnSplit::CopyToAll, OnSplit::Drop] {
                    for merge in [OnMerge::KeepIfAllEqual, OnMerge::Drop] {
                        for equal in [true, false] {
                            let pol = policy(modify, transform, split, merge);
                            let table = PolicyTable(BTreeMap::from([(AttributeKey(1), pol)]));
                            let inputs: AttributeMap = (1..=6)
                                .map(|n| {
                                    let v: &[u8] = if n == 6 && !equal { b"y" } else { b"x" };
                                    (id(n), vec![attr(1, v)])
                                })
                                .collect();
                            let mut h = History::new(
                                OperationId(1),
                                kind,
                                vec![],
                                vec![],
                                relations.clone(),
                                vec![],
                            );
                            let outcome = |from: u8, to: Vec<u8>, result| AttributeOutcome {
                                key: AttributeKey(1),
                                from: id(from),
                                to: to.into_iter().map(id).collect(),
                                result,
                            };
                            let (m_to, m_res) = match (kind, modify, transform) {
                                (OperationKind::Transform, _, OnTransform::Keep) => (vec![4], Kept),
                                (OperationKind::Transform, _, OnTransform::Drop) => {
                                    (vec![], Dropped)
                                }
                                (OperationKind::Transform, _, OnTransform::Recompute) => {
                                    (vec![4], Recomputed)
                                }
                                (_, OnModify::Keep, _) => (vec![4], Kept),
                                (_, OnModify::Drop, _) => (vec![], Dropped),
                            };
                            let (s_to, s_res) = match split {
                                OnSplit::CopyToAll => (vec![10, 11], Copied),
                                OnSplit::Drop => (vec![], Dropped),
                            };
                            let (g_to, g_res) = match (merge, equal) {
                                (OnMerge::Drop, _) => (vec![], Dropped),
                                (_, true) => (vec![12], Kept),
                                (_, false) => (vec![], Conflict),
                            };
                            h.attributes = vec![
                                outcome(1, s_to.clone(), s_res),
                                outcome(2, vec![2], Kept),
                                outcome(3, vec![], Dropped),
                                outcome(4, m_to.clone(), m_res),
                                outcome(5, g_to.clone(), g_res),
                                outcome(6, g_to.clone(), g_res),
                            ];
                            let mut outputs = AttributeMap::new();
                            for (from, to) in [(1, &s_to), (2, &vec![2]), (4, &m_to), (5, &g_to)] {
                                for t in to {
                                    let v: &[u8] = if m_res == Recomputed && from == 4 {
                                        b"new"
                                    } else {
                                        b"x"
                                    };
                                    outputs.entry(id(*t)).or_default().push(attr(1, v));
                                }
                            }
                            let issues = attributes::check(&inputs, &outputs, &h, &table);
                            assert_eq!(issues, vec![], "{kind:?} {pol:?} equal={equal}");
                            // Contradicting the policy is always reported.
                            let mut wrong = h.clone();
                            wrong.attributes[3].result = if m_res == Kept { Dropped } else { Kept };
                            assert!(
                                !attributes::check(&inputs, &outputs, &wrong, &table).is_empty()
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn attribute_issues_name_their_key_and_entity() {
    use attributes::AttributeIssueKind as A;
    use OutcomeResult::*;
    let keep = policy(
        OnModify::Keep,
        OnTransform::Keep,
        OnSplit::CopyToAll,
        OnMerge::KeepIfAllEqual,
    );
    let table = PolicyTable(BTreeMap::from([(AttributeKey(1), keep)]));
    let relations = vec![
        Relation::Modified {
            from: id(1),
            to: id(1),
        },
        Relation::Deleted { id: id(2) },
    ];
    let inputs: AttributeMap =
        BTreeMap::from([(id(1), vec![attr(1, b"x")]), (id(2), vec![attr(1, b"z")])]);
    let good = vec![
        AttributeOutcome {
            key: AttributeKey(1),
            from: id(1),
            to: vec![id(1)],
            result: Kept,
        },
        AttributeOutcome {
            key: AttributeKey(1),
            from: id(2),
            to: vec![],
            result: Dropped,
        },
    ];
    let h = |attributes: Vec<AttributeOutcome>| {
        History::new(
            OperationId(1),
            OperationKind::Extrude,
            vec![],
            vec![],
            relations.clone(),
            attributes,
        )
    };
    let outputs: AttributeMap = BTreeMap::from([(id(1), vec![attr(1, b"x")])]);
    let issues =
        |h: &History, inputs: &AttributeMap, outputs: &AttributeMap, table: &PolicyTable| {
            attributes::check(inputs, outputs, h, table)
                .iter()
                .map(|i| (i.kind, i.id))
                .collect::<Vec<_>>()
        };
    assert_eq!(issues(&h(good.clone()), &inputs, &outputs, &table), vec![]);
    assert_eq!(
        issues(&h(good.clone()), &inputs, &outputs, &PolicyTable::default()),
        // Without a policy nothing may carry the output attribute either.
        vec![
            (A::MissingPolicy, id(1)),
            (A::MissingPolicy, id(2)),
            (A::UntracedAttribute, id(1))
        ]
    );
    assert_eq!(
        issues(&h(good[..1].to_vec()), &inputs, &outputs, &table),
        vec![(A::MissingOutcome, id(2))]
    );
    let mut twice = good.clone();
    twice.push(good[1].clone());
    assert_eq!(
        issues(&h(twice), &inputs, &outputs, &table),
        vec![(A::DuplicateOutcome, id(2))]
    );
    let mut dangling = good.clone();
    dangling.push(AttributeOutcome {
        key: AttributeKey(1),
        from: id(3),
        to: vec![],
        result: Dropped,
    });
    assert_eq!(
        issues(&h(dangling), &inputs, &outputs, &table),
        vec![(A::DanglingOutcome, id(3))]
    );
    let changed: AttributeMap = BTreeMap::from([(id(1), vec![attr(1, b"q")])]);
    assert_eq!(
        issues(&h(good.clone()), &inputs, &changed, &table),
        vec![(A::ValueChanged, id(1))]
    );
    assert_eq!(
        issues(&h(good.clone()), &inputs, &AttributeMap::new(), &table),
        vec![(A::MissingAttribute, id(1))]
    );
    let extra: AttributeMap = BTreeMap::from([(id(1), vec![attr(1, b"x"), attr(2, b"n")])]);
    assert_eq!(
        issues(&h(good.clone()), &inputs, &extra, &table),
        vec![(A::UntracedAttribute, id(1))]
    );
    // An attribute left on an output after every outcome dropped it.
    let mut dropped = good.clone();
    dropped[0] = AttributeOutcome {
        key: AttributeKey(1),
        from: id(1),
        to: vec![],
        result: Dropped,
    };
    assert_eq!(
        issues(&h(dropped), &inputs, &outputs, &table),
        vec![
            (A::OutcomeContradictsPolicy, id(1)),
            (A::AttributeOnDropped, id(1))
        ]
    );
}
