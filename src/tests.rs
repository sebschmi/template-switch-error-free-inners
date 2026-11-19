use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use traitsequence::interface::Sequence;

use crate::MatchTable;

#[test]
fn reference_reference() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AGGGGAACCCCAA").unwrap();
    let query = VectorGenome::from_slice_u8(b"AAAAAAAA").unwrap();
    let matches = MatchTable::new(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        4,
    );

    for primary_index in 0..reference.len() - 3 {
        for secondary_rc_index in 0..reference.len() - 3 {
            assert_eq!(
                matches.has_reference_reference_match(primary_index, secondary_rc_index),
                (primary_index == 7 && secondary_rc_index == 8)
                    || (primary_index == 1 && secondary_rc_index == 2),
                "{primary_index}/{secondary_rc_index}"
            );
        }
    }

    for reference_kmer_index in 0..reference.len() - 3 {
        for query_kmer_index in 0..query.len() - 3 {
            assert!(!matches.has_reference_query_match(reference_kmer_index, query_kmer_index));
            assert!(!matches.has_query_reference_match(query_kmer_index, reference_kmer_index));
        }
    }

    for primary_index in 0..query.len() - 3 {
        for secondary_rc_index in 0..query.len() - 3 {
            assert!(!matches.has_query_query_match(primary_index, secondary_rc_index));
        }
    }
}

#[test]
fn reference_query() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AGGGGAA").unwrap();
    let query = VectorGenome::from_slice_u8(b"AACCCCA").unwrap();
    let matches = MatchTable::new(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        4,
    );

    for primary_index in 0..reference.len() - 3 {
        for secondary_rc_index in 0..reference.len() - 3 {
            assert!(!matches.has_reference_reference_match(primary_index, secondary_rc_index));
        }
    }

    for reference_kmer_index in 0..reference.len() - 3 {
        for query_kmer_index in 0..query.len() - 3 {
            assert_eq!(
                matches.has_reference_query_match(reference_kmer_index, query_kmer_index),
                reference_kmer_index == 1 && query_kmer_index == 1,
                "{reference_kmer_index}/{query_kmer_index}"
            );
            assert_eq!(
                matches.has_query_reference_match(query_kmer_index, reference_kmer_index),
                query_kmer_index == 2 && reference_kmer_index == 2,
                "{reference_kmer_index}/{query_kmer_index}"
            );
        }
    }

    for primary_index in 0..query.len() - 3 {
        for secondary_rc_index in 0..query.len() - 3 {
            assert!(!matches.has_query_query_match(primary_index, secondary_rc_index));
        }
    }
}

#[test]
fn sequences_shorter_than_minimum_length() {
    // Regression test: sequences shorter than `minimum_length` used to underflow the kmer count.
    for (reference, query) in [
        (&b"ACG"[..], &b"AACCCCAGGGGA"[..]),
        (&b"AACCCCAGGGGA"[..], &b"ACG"[..]),
        (&b"ACG"[..], &b"TT"[..]),
        (&b""[..], &b""[..]),
    ] {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(reference).unwrap();
        let query = VectorGenome::from_slice_u8(query).unwrap();
        let matches = MatchTable::new(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            4,
        );

        let reference_kmer_count = (reference.len() + 1).saturating_sub(4);
        let query_kmer_count = (query.len() + 1).saturating_sub(4);
        assert_eq!(matches.reference_kmer_count, reference_kmer_count);
        assert_eq!(matches.query_kmer_count, query_kmer_count);
        assert!(reference_kmer_count == 0 || query_kmer_count == 0);

        // The tables involving a short sequence are empty.
        assert!(matches.reference_query.is_empty());
        assert!(matches.query_reference.is_empty());
        if reference_kmer_count == 0 {
            assert!(matches.reference_reference.is_empty());
        }
        if query_kmer_count == 0 {
            assert!(matches.query_query.is_empty());
        }

        // The table of the long sequence is still computed correctly.
        for primary_index in 0..reference_kmer_count {
            for secondary_rc_index in 0..reference_kmer_count {
                assert_eq!(
                    matches.has_reference_reference_match(primary_index, secondary_rc_index),
                    (primary_index == 7 && secondary_rc_index == 6)
                        || (primary_index == 2 && secondary_rc_index == 1),
                    "{primary_index}/{secondary_rc_index}"
                );
            }
        }
        for primary_index in 0..query_kmer_count {
            for secondary_rc_index in 0..query_kmer_count {
                assert_eq!(
                    matches.has_query_query_match(primary_index, secondary_rc_index),
                    (primary_index == 7 && secondary_rc_index == 6)
                        || (primary_index == 2 && secondary_rc_index == 1),
                    "{primary_index}/{secondary_rc_index}"
                );
            }
        }
    }
}

#[test]
fn query_query() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AAAAAAA").unwrap();
    let query = VectorGenome::from_slice_u8(b"AACCCCAGGGGA").unwrap();
    let matches = MatchTable::new(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        4,
    );

    for primary_index in 0..reference.len() - 3 {
        for secondary_rc_index in 0..reference.len() - 3 {
            assert!(!matches.has_reference_reference_match(primary_index, secondary_rc_index));
        }
    }

    for reference_kmer_index in 0..reference.len() - 3 {
        for query_kmer_index in 0..query.len() - 3 {
            assert!(!matches.has_reference_query_match(reference_kmer_index, query_kmer_index));
            assert!(!matches.has_query_reference_match(query_kmer_index, reference_kmer_index));
        }
    }

    for primary_index in 0..query.len() - 3 {
        for secondary_rc_index in 0..query.len() - 3 {
            assert_eq!(
                matches.has_query_query_match(primary_index, secondary_rc_index),
                (primary_index == 7 && secondary_rc_index == 6)
                    || (primary_index == 2 && secondary_rc_index == 1),
                "{primary_index}/{secondary_rc_index}"
            );
        }
    }
}
