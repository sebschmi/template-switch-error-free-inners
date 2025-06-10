//! Compute all error-free template switch inners for a pair of genome strings.

#![warn(missing_docs)]

use bitvec::vec::BitVec;
use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use log::debug;
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

        debug!("Converting genomes to strings");
        let reference_rc =
            VectorGenome::<AlphabetType>::from_iter(reference.reverse_complement_iter())
                .as_string();
        let query_rc =
            VectorGenome::<AlphabetType>::from_iter(query.reverse_complement_iter()).as_string();
        let reference = reference.as_string();
        let query = query.as_string();

        debug!("Computing indexes");
        let reference_index = SuffixTable::new(&reference);
        let query_index = SuffixTable::new(&query);

        debug!("Initialising bitvectors");
        let reference_kmer_count = reference.len() - minimum_length + 1;
        let query_kmer_count = query.len() - minimum_length + 1;

        let reference_reference =
            BitVec::repeat(false, reference_kmer_count * reference_kmer_count);
        let reference_query = BitVec::repeat(false, reference_kmer_count * query_kmer_count);
        let query_reference = BitVec::repeat(false, query_kmer_count * reference_kmer_count);
        let query_query = BitVec::repeat(false, query_kmer_count * query_kmer_count);

        debug!("Finding matches");
        let reference_rc_character_offsets: Vec<_> = reference_rc
            .char_indices()
            .map(|(index, _)| index)
            .collect();
        let query_rc_character_offsets: Vec<_> =
            query_rc.char_indices().map(|(index, _)| index).collect();

        for reference_rc_kmer_index in 0..reference_kmer_count {
            todo!();
        }

        for query_rc_kmer_index in 0..query_kmer_count {
            todo!();
        }

        Self {
            reference_reference,
            reference_query,
            query_reference,
            query_query,
        }
    }
}
