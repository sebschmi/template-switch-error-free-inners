//! Compute all error-free template switch inners for a pair of genome strings.

#![warn(missing_docs)]

use std::iter;

use bitvec::vec::BitVec;
use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use log::debug;
use suffix::SuffixTable;

#[cfg(test)]
mod tests;

/// A table of all error-free template switch inner entry points for a pair of genome strings.
pub struct MatchTable {
    reference_reference: BitVec,
    reference_query: BitVec,
    query_reference: BitVec,
    query_query: BitVec,
    reference_kmer_count: usize,
    query_kmer_count: usize,
}

impl MatchTable {
    /// Compute all error-free template switch inner entry points for a pair of genome strings.
    ///
    /// The inners must have the given minimum length.
    pub fn new<
        AlphabetType: Alphabet,
        GenomeSubsequence: GenomeSequence<AlphabetType, GenomeSubsequence> + ?Sized,
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

        let mut reference_reference =
            BitVec::repeat(false, reference_kmer_count * reference_kmer_count);
        let mut reference_query = BitVec::repeat(false, reference_kmer_count * query_kmer_count);
        let mut query_reference = BitVec::repeat(false, query_kmer_count * reference_kmer_count);
        let mut query_query = BitVec::repeat(false, query_kmer_count * query_kmer_count);

        debug!("Finding matches");
        let reference_rc_character_offsets: Vec<_> = reference_rc
            .char_indices()
            .map(|(index, _)| index)
            .chain(iter::once(reference_rc.len()))
            .collect();
        let query_rc_character_offsets: Vec<_> = query_rc
            .char_indices()
            .map(|(index, _)| index)
            .chain(iter::once(query_rc.len()))
            .collect();

        for reference_rc_kmer_index in 0..reference_kmer_count {
            let reference_rc_kmer = &reference_rc[reference_rc_character_offsets
                [reference_rc_kmer_index]
                ..reference_rc_character_offsets[reference_rc_kmer_index + minimum_length]];

            for reference_kmer_index in reference_index.positions(reference_rc_kmer) {
                let reference_kmer_index = usize::try_from(*reference_kmer_index).unwrap();
                reference_reference.set(
                    reference_kmer_index * reference_kmer_count + reference_rc_kmer_index,
                    true,
                );
            }

            for query_kmer_index in query_index.positions(reference_rc_kmer) {
                let query_kmer_index = usize::try_from(*query_kmer_index).unwrap();
                query_reference.set(
                    query_kmer_index * reference_kmer_count + reference_rc_kmer_index,
                    true,
                );
            }
        }

        for query_rc_kmer_index in 0..query_kmer_count {
            let query_rc_kmer = &query_rc[query_rc_character_offsets[query_rc_kmer_index]
                ..query_rc_character_offsets[query_rc_kmer_index + minimum_length]];

            for reference_kmer_index in reference_index.positions(query_rc_kmer) {
                let reference_kmer_index = usize::try_from(*reference_kmer_index).unwrap();
                reference_query.set(
                    reference_kmer_index * query_kmer_count + query_rc_kmer_index,
                    true,
                );
            }

            for query_kmer_index in query_index.positions(query_rc_kmer) {
                let query_kmer_index = usize::try_from(*query_kmer_index).unwrap();
                query_query.set(
                    query_kmer_index * query_kmer_count + query_rc_kmer_index,
                    true,
                );
            }
        }

        Self {
            reference_reference,
            reference_query,
            query_reference,
            query_query,
            reference_kmer_count,
            query_kmer_count,
        }
    }

    /// Returns `true` if the reference kmer at `primary_index` matches the kmer in the reverse-complemented reference at `secondary_rc_index`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
    /// use compact_genome::implementation::vec_sequence::VectorGenome;
    /// use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    /// use template_switch_error_free_inners::MatchTable;
    ///
    /// let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AGGGGAACCCCAA").unwrap();
    /// let query = VectorGenome::from_slice_u8(b"AAAAAAAA").unwrap();
    /// let matches = MatchTable::new(
    ///     reference.as_genome_subsequence(),
    ///     query.as_genome_subsequence(),
    ///     4,
    /// );
    ///
    /// assert!(matches.has_reference_reference_match(1, 2));
    /// assert!(matches.has_reference_reference_match(7, 8));
    /// ```
    pub fn has_reference_reference_match(
        &self,
        primary_index: usize,
        secondary_rc_index: usize,
    ) -> bool {
        debug_assert!(primary_index < self.reference_kmer_count);
        debug_assert!(secondary_rc_index < self.reference_kmer_count);
        self.reference_reference[primary_index * self.reference_kmer_count + secondary_rc_index]
    }

    /// Returns `true` if the reference kmer at `primary_index` matches the kmer in the reverse-complemented query at `secondary_rc_index`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
    /// use compact_genome::implementation::vec_sequence::VectorGenome;
    /// use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    /// use template_switch_error_free_inners::MatchTable;
    ///
    /// let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AGGGGAAA").unwrap();
    /// let query = VectorGenome::from_slice_u8(b"ACCCCAA").unwrap();
    /// let matches = MatchTable::new(
    ///     reference.as_genome_subsequence(),
    ///     query.as_genome_subsequence(),
    ///     4,
    /// );
    ///
    /// assert!(matches.has_reference_query_match(1, 2));
    /// ```
    pub fn has_reference_query_match(
        &self,
        primary_index: usize,
        secondary_rc_index: usize,
    ) -> bool {
        debug_assert!(primary_index < self.reference_kmer_count);
        debug_assert!(secondary_rc_index < self.query_kmer_count);
        self.reference_query[primary_index * self.query_kmer_count + secondary_rc_index]
    }

    /// Returns `true` if the query kmer at `primary_index` matches the kmer in the reverse-complemented reference at `secondary_rc_index`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
    /// use compact_genome::implementation::vec_sequence::VectorGenome;
    /// use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    /// use template_switch_error_free_inners::MatchTable;
    ///
    /// let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AGGGGAAA").unwrap();
    /// let query = VectorGenome::from_slice_u8(b"ACCCCAA").unwrap();
    /// let matches = MatchTable::new(
    ///     reference.as_genome_subsequence(),
    ///     query.as_genome_subsequence(),
    ///     4,
    /// );
    ///
    /// assert!(matches.has_query_reference_match(1, 3));
    /// ```
    pub fn has_query_reference_match(
        &self,
        primary_index: usize,
        secondary_rc_index: usize,
    ) -> bool {
        debug_assert!(primary_index < self.query_kmer_count);
        debug_assert!(secondary_rc_index < self.reference_kmer_count);
        self.query_reference[primary_index * self.reference_kmer_count + secondary_rc_index]
    }

    /// Returns `true` if the query kmer at `primary_index` matches the kmer in the reverse-complemented query at `secondary_rc_index`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
    /// use compact_genome::implementation::vec_sequence::VectorGenome;
    /// use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    /// use template_switch_error_free_inners::MatchTable;
    ///
    /// let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"AAAAAA").unwrap();
    /// let query = VectorGenome::from_slice_u8(b"AAGGGGACCCCA").unwrap();
    /// let matches = MatchTable::new(
    ///     reference.as_genome_subsequence(),
    ///     query.as_genome_subsequence(),
    ///     4,
    /// );
    ///
    /// assert!(matches.has_query_query_match(2, 1));
    /// assert!(matches.has_query_query_match(7, 6));
    /// ```
    pub fn has_query_query_match(&self, primary_index: usize, secondary_rc_index: usize) -> bool {
        debug_assert!(primary_index < self.query_kmer_count);
        debug_assert!(secondary_rc_index < self.query_kmer_count);
        self.query_query[primary_index * self.query_kmer_count + secondary_rc_index]
    }
}
