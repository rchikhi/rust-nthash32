// 32 bit hashes version, used to compare against mine in rust-seq2kminmers as well as nthash.hpp
//    of the avx port
//
//! ntHash is a hash function tuned for genomic data.
//! It performs best when calculating hash values for adjacent k-mers in
//! an input sequence, operating an order of magnitude faster than the best
//! performing alternatives in typical use cases.
//!
//! [Scientific article with more details](https://doi.org/10.1093/bioinformatics/btw397)
//!
//! [Original implementation in C++](https://github.com/bcgsc/ntHash/)
//!
//! This crate is based on ntHash [1.0.4](https://github.com/bcgsc/ntHash/releases/tag/v1.0.4).
//!

mod error;

pub use crate::error::{Error, Result};

pub(crate) const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;

const SHIFT :u64 =0;

const H_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474 >> SHIFT;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c >> SHIFT;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324 >> SHIFT;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456 >> SHIFT;
    lookup[b'N' as usize] = 0;
    lookup
};

const RC_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456 >> SHIFT;
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324 >> SHIFT;
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c >> SHIFT;
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474 >> SHIFT;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> u32 {
    let val = H_LOOKUP[c as usize] as u32;
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }
    val
}

#[inline(always)]
fn rc(nt: u8) -> u32 {
    let val = RC_LOOKUP[nt as usize] as u32;
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", nt as char)
    }
    val
}

/// An efficient iterator for calculating hashes for genomic sequences.
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<u32>`.
/// ```
///     # use nthash32::Result;
///     use nthash32::NtHashIterator;
///
///     # fn main() -> Result<()> {
///     let seq = b"ACTGC";
///     let iter = NtHashIterator::new(seq, 3)?;
///     let hashes: Vec<u32> = iter.collect();
///     assert_eq!(hashes,
///                vec![0x9b1eda9a185413ce, 0x9f6acfa2235b86fc, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
/// or, in one line:
/// ```
///     # use nthash32::Result;
///     use nthash32::NtHashIterator;
///
///     # fn main() -> Result<()> {
///     assert_eq!(NtHashIterator::new(b"ACTGC", 3)?.collect::<Vec<u32>>(),
///                vec![0x9b1eda9a185413ce, 0x9f6acfa2235b86fc, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
#[derive(Debug)]
pub struct NtHashIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u32,
    rh: u32,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashIterator<'a> {
    /// Creates a new NtHashIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize) -> Result<NtHashIterator<'a>> {
        if k > seq.len() {
            return Err(Error::KSizeOutOfRange {
                ksize: k,
                seq_size: seq.len(),
            });
        }
        if k > MAXIMUM_K_SIZE {
            return Err(Error::KSizeTooBig(k));
        }
        let mut fh = 0;
        for (i, v) in seq[0..k].iter().enumerate() {
            fh ^= h(*v).rotate_left((k - i - 1) as u32);
        }

        let mut rh = 0;
        for (i, v) in seq[0..k].iter().rev().enumerate() {
            rh ^= rc(*v).rotate_left((k - i - 1) as u32);
        }

        Ok(NtHashIterator {
            seq,
            k,
            fh,
            rh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        if self.current_idx == self.max_idx {
            return None;
        };

        if self.current_idx != 0 {
            let i = self.current_idx - 1;
            let seqi = self.seq[i];
            let seqk = self.seq[i + self.k];

            self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);

            self.rh = self.rh.rotate_right(1)
                ^ rc(seqi).rotate_right(1)
                ^ rc(seqk).rotate_left(self.k as u32 - 1);
        }
        
        //println!(" nthash-luiz i {} fh {:#x?} rh {:#x?}",self.current_idx, self.fh,self.rh);

        self.current_idx += 1;
        Some(u32::min(self.rh, self.fh))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for NtHashIterator<'a> {}

