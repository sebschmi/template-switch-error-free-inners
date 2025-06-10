//! Compute all error-free template switch inners for a pair of genome strings.

#![warn(missing_docs)]

use bitvec::vec::BitVec;
use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use suffix::SuffixTable;

/// A table of all error-free template switch inner entry points for a pair of genome strings.
pub struct MatchTable {
    reference_reference: BitVec,
    reference_query: BitVec,
    query_reference: BitVec,
    query_query: BitVec,
}

impl MatchTable {
    /// Compute all error-free template switch inner entry points for a pair of genome strings.
    ///
    /// The inners must have the given minimum length.
    pub fn new<
        AlphabetType: Alphabet,
        GenomeSubsequence: GenomeSequence<AlphabetType, GenomeSubsequence>,
    >(
        reference: &GenomeSubsequence,
        query: &GenomeSubsequence,
        minimum_length: usize,
    ) -> Self {
        assert!(minimum_length > 0);

        let reference_rc =
            VectorGenome::<AlphabetType>::from_iter(reference.reverse_complement_iter())
                .as_string();
        let query_rc =
            VectorGenome::<AlphabetType>::from_iter(query.reverse_complement_iter()).as_string();
        let reference = reference.as_string();
        let query = query.as_string();

        todo!()
    }
}
