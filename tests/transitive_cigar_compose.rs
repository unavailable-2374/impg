//! Integration tests for CIGAR-precise hop≥1 transitive BFS (the
//! `impg depth --use-BFS` core). Validates that with `store_cigar = true`
//! (carry mode), every transitive hit reports the ANCHOR as its `overlap.2`
//! target with a synthesized anchor→query CIGAR, and that the synthesized
//! CIGAR matches a hand-composed expectation through a 2-hop chain.
//!
//! See `notes/PLAN_hop_ge1_cigar_precise_IMPL.md`.

use impg::alignment_record::AlignmentRecord;
use impg::impg::{CigarOp, Impg};
use impg::seqidx::SequenceIndex;
use std::fs::File;
use std::io::Write;
use std::num::NonZeroUsize;
use tempfile::TempDir;

/// Build a single-file Impg from PAF lines written to a temp dir.
fn build_impg(lines: &[&str]) -> Impg {
    let tmp = TempDir::new().expect("tempdir");
    let paf_path = tmp.path().join("aln.paf");
    {
        let mut f = File::create(&paf_path).expect("create paf");
        for l in lines {
            writeln!(f, "{}", l).expect("write");
        }
    }
    let mut seq_index = SequenceIndex::new();
    let file = File::open(&paf_path).expect("open paf");
    let records: Vec<AlignmentRecord> = impg::paf::parse_paf_file(
        paf_path.to_str().unwrap(),
        file,
        NonZeroUsize::new(1).unwrap(),
        &mut seq_index,
    )
    .expect("parse paf");
    let records_by_file = vec![(records, paf_path.to_string_lossy().into_owned())];
    // Keep the TempDir alive for the lifetime of the index by leaking it: the
    // CIGARs are read lazily from the PAF file on disk during projection.
    std::mem::forget(tmp);
    Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true)
        .expect("from_multi_alignment_records")
}

fn cigar_string(ops: &[CigarOp]) -> String {
    // Merge adjacent equal ops (new_run chunking) into a readable string.
    let mut out = String::new();
    let mut run_len: i64 = 0;
    let mut run_op = '\0';
    for o in ops {
        if o.op() == run_op {
            run_len += o.len() as i64;
        } else {
            if run_op != '\0' {
                out.push_str(&format!("{}{}", run_len, run_op));
            }
            run_op = o.op();
            run_len = o.len() as i64;
        }
    }
    if run_op != '\0' {
        out.push_str(&format!("{}{}", run_len, run_op));
    }
    out
}

/// Walk one alignment interval (target forward from `t_start`, query oriented by
/// `q_first`/`q_last`) and return a map: target_pos -> aligned query base. This
/// is the trusted base-level projection used to build a 2-hop ground truth,
/// independent of the transitive composition under test.
fn base_map(
    t_start: i64,
    q_first: i64,
    q_last: i64,
    ops: &[CigarOp],
) -> std::collections::HashMap<i64, i64> {
    use std::collections::HashMap;
    let fwd = q_first <= q_last;
    let mut t = t_start;
    let mut q = if fwd { q_first } else { q_last };
    let mut m = HashMap::new();
    for o in ops {
        let n = o.len() as i64;
        match o.op() {
            '=' | 'X' | 'M' => {
                for _ in 0..n {
                    let qb = if fwd { q } else { q - 1 };
                    m.insert(t, qb);
                    t += 1;
                    q += if fwd { 1 } else { -1 };
                }
            }
            'D' => t += n,
            'I' => q += if fwd { n } else { -n },
            _ => {}
        }
    }
    m
}

#[test]
fn forward_two_hop_compose_targets_anchor() {
    // A↔B identity (100=), B↔C has a 10bp deletion in C at B[40,50).
    // Expect: A[0,100) projects to C[0,90) with CIGAR "40=10D50=", target=A.
    let impg = build_impg(&[
        // query=B target=A, identity
        "B\t100\t0\t100\t+\tA\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        // query=C target=B, 10bp deletion in C (D consumes target B only)
        "C\t100\t0\t90\t+\tB\t100\t0\t100\t90\t90\t60\tcg:Z:40=10D50=",
    ]);
    let a = impg.seq_index.get_id("A").expect("A id");
    let c = impg.seq_index.get_id("C").expect("C id");

    let results = impg.query_transitive_bfs(
        a, 0, 100, None, 0, 1, 0, None, /*store_cigar=*/ true, None, None, false, None,
    );

    // Find the hop to C.
    let to_c: Vec<_> = results.iter().filter(|(q, _, _)| q.metadata == c).collect();
    assert!(!to_c.is_empty(), "no transitive hit reached C");

    for (q, cig, t) in &to_c {
        // overlap.2 (target) must be the ANCHOR A, not the hub B (§5.4).
        assert_eq!(
            t.metadata, a,
            "carry-mode transitive hit must report anchor A as target, got {}",
            t.metadata
        );
        // Anchor coords forward and within the queried region.
        assert!(t.first < t.last && t.first >= 0 && t.last <= 100);
        // The synthesized CIGAR is the hand-composed A→C.
        assert_eq!(cigar_string(cig), "40=10D50=");
        // C span: forward [0,90).
        let (c_lo, c_hi) = (q.first.min(q.last), q.first.max(q.last));
        assert_eq!((c_lo, c_hi), (0, 90));
        // Target (anchor) span covers all of A.
        assert_eq!((t.first, t.last), (0, 100));
    }
}

#[test]
fn non_carry_mode_targets_hub_not_anchor() {
    // With store_cigar=false the legacy path must be unchanged: the hop-2 hit
    // reports the immediate hub B as its target, not the anchor A.
    let impg = build_impg(&[
        "B\t100\t0\t100\t+\tA\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        "C\t100\t0\t90\t+\tB\t100\t0\t100\t90\t90\t60\tcg:Z:40=10D50=",
    ]);
    let a = impg.seq_index.get_id("A").unwrap();
    let b = impg.seq_index.get_id("B").unwrap();
    let c = impg.seq_index.get_id("C").unwrap();

    let results = impg.query_transitive_bfs(
        a, 0, 100, None, 0, 1, 0, None, /*store_cigar=*/ false, None, None, false, None,
    );
    let to_c: Vec<_> = results.iter().filter(|(q, _, _)| q.metadata == c).collect();
    assert!(!to_c.is_empty(), "no transitive hit reached C");
    for (_, _, t) in &to_c {
        assert_eq!(
            t.metadata, b,
            "legacy path must report hub B as target for the hop-2 hit"
        );
    }
}

#[test]
fn dfs_two_hop_compose_matches_bfs() {
    // The DFS carry path must synthesize the same anchor→C CIGAR as BFS.
    let impg = build_impg(&[
        "B\t100\t0\t100\t+\tA\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        "C\t100\t0\t90\t+\tB\t100\t0\t100\t90\t90\t60\tcg:Z:40=10D50=",
    ]);
    let a = impg.seq_index.get_id("A").unwrap();
    let c = impg.seq_index.get_id("C").unwrap();

    let results = impg.query_transitive_dfs(
        a, 0, 100, None, 0, 1, 0, None, true, None, None, false, None,
    );
    let to_c: Vec<_> = results.iter().filter(|(q, _, _)| q.metadata == c).collect();
    assert!(!to_c.is_empty(), "no transitive hit reached C via DFS");
    for (q, cig, t) in &to_c {
        assert_eq!(t.metadata, a, "DFS carry hit must target the anchor");
        assert_eq!(cigar_string(cig), "40=10D50=");
        assert_eq!((q.first.min(q.last), q.first.max(q.last)), (0, 90));
        assert_eq!((t.first, t.last), (0, 100));
    }
}

#[test]
fn reverse_hop_compose_directions() {
    // A↔B identity forward. B↔C reverse: C aligns to B on the minus strand.
    // The composite A→C must therefore be reverse (first > last on the C side).
    let impg = build_impg(&[
        "B\t100\t0\t100\t+\tA\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        // query=C target=B, reverse strand, identity 100=
        "C\t100\t0\t100\t-\tB\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
    ]);
    let a = impg.seq_index.get_id("A").unwrap();
    let c = impg.seq_index.get_id("C").unwrap();

    let results = impg.query_transitive_bfs(
        a, 0, 100, None, 0, 1, 0, None, true, None, None, false, None,
    );
    let to_c: Vec<_> = results.iter().filter(|(q, _, _)| q.metadata == c).collect();
    assert!(!to_c.is_empty(), "no transitive hit reached C");
    for (q, _cig, t) in &to_c {
        // Anchor target always forward.
        assert_eq!(t.metadata, a);
        assert!(t.first < t.last);
        // C↔A is reverse, so the query (C) interval must be stored reversed
        // (first > last), the convention depth's CigarEntry relies on (§5.4).
        assert!(
            q.first > q.last,
            "reverse composite must store C interval reversed: got first={}, last={}",
            q.first,
            q.last
        );
    }
}

#[test]
fn reverse_hop_with_indel_is_single_base_exact() {
    // Regression for the reverse-strand + indel composition bug: A↔B forward
    // identity, B↔C REVERSE with a 2bp deletion in C. The composed A→C must map
    // every anchor base to the exact same C base as chaining the two trusted
    // single-hop projections (position-chain ground truth). Before the fix the
    // indel was misplaced because the carried CIGAR's B axis was reversed
    // relative to the pairwise B axis.
    let impg = build_impg(&[
        "B\t20\t0\t20\t+\tA\t20\t0\t20\t20\t20\t60\tcg:Z:20=",
        // query=C target=B, reverse, 2bp deletion in C: target B=20, query C=18
        "C\t18\t0\t18\t-\tB\t20\t0\t20\t18\t18\t60\tcg:Z:8=2D10=",
    ]);
    let a = impg.seq_index.get_id("A").unwrap();
    let b = impg.seq_index.get_id("B").unwrap();
    let c = impg.seq_index.get_id("C").unwrap();

    let res_a =
        impg.query_transitive_bfs(a, 0, 20, None, 0, 1, 0, None, true, None, None, false, None);
    let res_b =
        impg.query_transitive_bfs(b, 0, 20, None, 0, 1, 0, None, true, None, None, false, None);

    // Trusted single hops: A→B (from anchor A) and B→C (from anchor B).
    let (ts_ab, q_ab, ops_ab) = res_a
        .iter()
        .find(|(q, _, t)| q.metadata == b && t.metadata == a)
        .map(|(q, cg, t)| (t.first.min(t.last), (q.first, q.last), cg.clone()))
        .expect("A→B hop");
    let (ts_bc, q_bc, ops_bc) = res_b
        .iter()
        .find(|(q, _, t)| q.metadata == c && t.metadata == b)
        .map(|(q, cg, t)| (t.first.min(t.last), (q.first, q.last), cg.clone()))
        .expect("B→C hop");

    // Ground truth: A-base -> B-base -> C-base.
    let m_ab = base_map(ts_ab, q_ab.0, q_ab.1, &ops_ab);
    let m_bc = base_map(ts_bc, q_bc.0, q_bc.1, &ops_bc);
    let mut gt = std::collections::HashMap::new();
    for (&apos, &bpos) in &m_ab {
        if let Some(&cpos) = m_bc.get(&bpos) {
            gt.insert(apos, cpos);
        }
    }

    // Device under test: composed A→C.
    let (ts_ac, q_ac, ops_ac) = res_a
        .iter()
        .find(|(q, _, t)| q.metadata == c && t.metadata == a)
        .map(|(q, cg, t)| (t.first.min(t.last), (q.first, q.last), cg.clone()))
        .expect("A→C composed hop");
    let mt = base_map(ts_ac, q_ac.0, q_ac.1, &ops_ac);

    assert!(!gt.is_empty() && !mt.is_empty());
    // A→C is two strand flips from... actually F then R → reverse; C interval reversed.
    assert!(
        q_ac.0 > q_ac.1,
        "A→C through one reverse hop must be reverse"
    );
    let mut compared = 0;
    for (&apos, &cpos) in &gt {
        if let Some(&mtc) = mt.get(&apos) {
            assert_eq!(
                mtc, cpos,
                "single-base mismatch at anchor {}: composed={} groundtruth={}",
                apos, mtc, cpos
            );
            compared += 1;
        }
    }
    assert_eq!(
        compared,
        gt.len(),
        "every ground-truth anchor base must be present in the composed map"
    );
}

#[test]
fn mixed_strand_three_hop_strand_accumulates() {
    // A→B forward, B→C reverse, C→D reverse. Composite A→D strand =
    // F XOR R XOR R = Forward, so the two reverses must cancel and the final D
    // interval must be stored FORWARD (first < last). This exercises the strand
    // XOR accumulation across two transposes — risk #1 in the plan (§8.1).
    let impg = build_impg(&[
        "B\t100\t0\t100\t+\tA\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        "C\t100\t0\t100\t-\tB\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
        "D\t100\t0\t100\t-\tC\t100\t0\t100\t100\t100\t60\tcg:Z:100=",
    ]);
    let a = impg.seq_index.get_id("A").unwrap();
    let c = impg.seq_index.get_id("C").unwrap();
    let d = impg.seq_index.get_id("D").unwrap();

    let results = impg.query_transitive_bfs(
        a, 0, 100, None, 0, 1, 0, None, true, None, None, false, None,
    );

    // C is one reverse hop from A → C interval must be reversed.
    for (q, _, t) in results.iter().filter(|(q, _, _)| q.metadata == c) {
        assert_eq!(t.metadata, a);
        assert!(q.first > q.last, "A→C (single reverse) must be reversed");
    }
    // D is two reverse hops from A → strand cancels → forward.
    let to_d: Vec<_> = results.iter().filter(|(q, _, _)| q.metadata == d).collect();
    assert!(!to_d.is_empty(), "no transitive hit reached D");
    for (q, _, t) in &to_d {
        assert_eq!(t.metadata, a, "3-hop hit must still target the anchor");
        assert!(t.first < t.last);
        assert!(
            q.first < q.last,
            "A→D (two reverses cancel) must be forward: got first={}, last={}",
            q.first,
            q.last
        );
        assert_eq!((q.first.min(q.last), q.first.max(q.last)), (0, 100));
        assert_eq!((t.first, t.last), (0, 100));
    }
}
