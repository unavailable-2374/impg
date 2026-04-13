//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!
//! The data structure used a fairly standard interval tree, but with nodes stored
//! in van Emde Boas layout, which improves average cache locality, and thus
//! query performance. The downside is that building the tree is more expensive
//! so a relatively large number of queries needs to made for it to break even.
//!
//! The data structure `COITree` is constructed with an array of `IntervalNode`
//! structs which store integer, end-inclusive intervals along with associated
//! metadata. The tree can be queried directly for coverage or overlaps, or
//! through the intermediary `SortedQuerenty` which keeps track of some state
//! to accelerate overlaping queries.

mod interval;
pub use interval::*;

mod nosimd;
pub use nosimd::*;

// nosimd feature: skip all SIMD implementations entirely
#[cfg(not(feature = "nosimd"))]
#[cfg(all(target_feature = "avx2"))]
mod avx;
#[cfg(not(feature = "nosimd"))]
#[cfg(all(target_feature = "avx2"))]
pub use avx::*;

#[cfg(not(feature = "nosimd"))]
#[cfg(all(target_feature = "neon"))]
mod neon;
#[cfg(not(feature = "nosimd"))]
#[cfg(all(target_feature = "neon"))]
pub use neon::*;

// When nosimd feature is active, always alias COITree to BasicCOITree
#[cfg(feature = "nosimd")]
pub type COITree<T, I> = BasicCOITree<T, I>;
#[cfg(feature = "nosimd")]
pub type COITreeSortedQuerent<'a, T, I> = BasicSortedQuerent<'a, T, I>;

// When nosimd feature is NOT active, pick based on target features
#[cfg(all(not(feature = "nosimd"), not(target_feature = "avx2"), not(target_feature = "neon")))]
pub type COITree<T, I> = BasicCOITree<T, I>;
#[cfg(all(not(feature = "nosimd"), not(target_feature = "avx2"), not(target_feature = "neon")))]
pub type COITreeSortedQuerent<'a, T, I> = BasicSortedQuerent<'a, T, I>;

#[cfg(all(not(feature = "nosimd"), target_feature = "avx2"))]
pub type COITree<T, I> = AVXCOITree<T, I>;
#[cfg(all(not(feature = "nosimd"), target_feature = "avx2"))]
pub type COITreeSortedQuerent<'a, T, I> = AVXSortedQuerent<'a, T, I>;

#[cfg(all(not(feature = "nosimd"), target_feature = "neon"))]
pub type COITree<T, I> = NeonCOITree<T, I>;
#[cfg(all(not(feature = "nosimd"), target_feature = "neon"))]
pub type COITreeSortedQuerent<'a, T, I> = NeonSortedQuerent<'a, T, I>;
