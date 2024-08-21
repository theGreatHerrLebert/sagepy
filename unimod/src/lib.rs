pub mod unimod {
    pub mod modification_atomic_composition;
    pub mod title_to_unimod_id;
    pub mod unimod_quantized;
    pub mod unimod_to_mass;

    // Re-exporting functions to the parent module for easier access
    pub use modification_atomic_composition::modification_atomic_composition;
    pub use title_to_unimod_id::title_to_unimod_id;
    pub use unimod_quantized::{quanzie_mass, quantized_mass_to_unimod};
    pub use unimod_to_mass::{unimod_modifications_mass, unimod_modifications_mass_numerical};
}