#[macro_use]
extern crate quickcheck;

use quickcheck::{Arbitrary, Gen};
use rand::Rng;

use nthash::{nthash, NtHashIterator};

#[test]
fn iter_cmp() {
    let ksize = 5;
    for s in &vec!["TGCAG", "ACGTC", "ACGTCGTCAGTCGATGCAGT", "ACGTCGANNGTA"] {
        let seq = s.as_bytes();
        let iter = NtHashIterator::new(seq, ksize).unwrap();
        println!("{:?}", s);
        assert_eq!(nthash(seq, ksize), iter.collect::<Vec<u64>>());
    }
}

#[should_panic(expected = "Non-ACGTN nucleotide encountered! E")]
#[test]
fn panic_non_acgtn() {
    let ksize: usize = 2;
    let sequences = "TGCAGNE";
    let iter = NtHashIterator::new(sequences.as_bytes(), ksize).unwrap();
    let _: Vec<u64> = iter.collect();
}

#[test]
fn out_of_range_ksize_wont_panic() {
    let ksize: usize = 10;
    let sequences = "TGCAG";
    let err = NtHashIterator::new(sequences.as_bytes(), ksize).unwrap_err();
    assert_eq!(
        err.to_string(),
        "K size 10 is out of range for the given sequence size 5"
    );
}

#[cfg(target_pointer_width = "64")]
#[test]
#[ignore]
fn big_ksize_wont_panic() {
    let ksize: usize = (u64::from(u32::max_value()) + 1) as usize;
    let repetitions: usize = ((f64::from(u32::max_value()) + 1.0) / 5.0).ceil() as usize;
    let sequences = "TGCAG".repeat(repetitions);
    let err = NtHashIterator::new(sequences.as_bytes(), ksize).unwrap_err();
    assert_eq!(
        err.to_string(),
        "K size 4294967296 cannot exceed the size of a u32 4294967295"
    );
}

#[derive(Clone, Debug)]
struct Seq(String);

impl Arbitrary for Seq {
    fn arbitrary<G: Gen>(g: &mut G) -> Seq {
        let choices = ['A', 'C', 'G', 'T', 'N'];
        let size = {
            let s = g.size();
            g.gen_range(0, s)
        };
        let mut s = String::with_capacity(size);
        for _ in 0..size {
            s.push(*g.choose(&choices).expect("Not a valid nucleotide"));
        }
        Seq { 0: s }
    }
}

quickcheck! {
  fn oracle_quickcheck(s: Seq) -> bool {
     let seq = s.0.as_bytes();
     (1..(seq.len())).all(|ksize| {
       let iter = NtHashIterator::new(seq, ksize).unwrap();
       nthash(seq, ksize) == iter.collect::<Vec<u64>>()
     })
  }
}
