#![allow(clippy::op_ref)]

pub mod accumulation;
pub mod backend;
pub mod frontend;
pub mod pcs;
pub mod piop;
pub mod poly;
pub mod util;
pub mod logging;

pub use halo2_curves;

pub use util::{
    expression::{CommonPolynomial, Expression, Query, Rotation},
    Deserialize, DeserializeOwned, Serialize,
};

#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    InvalidSumcheck(String),
    InvalidPcsParam(String),
    InvalidPcsOpen(String),
    InvalidSnark(String),
    Serialization(String),
    Transcript(std::io::ErrorKind, String),
}
