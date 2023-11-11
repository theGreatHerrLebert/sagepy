use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::py_mass::PyTolerance;
use sage_core::spectrum::{
    Deisotoped, Peak, Precursor, ProcessedSpectrum, RawSpectrum, Representation, SpectrumProcessor,
};

#[pyclass]
#[derive(Clone)]
pub struct PyRepresentation {
    pub inner: Representation,
}

#[pymethods]
impl PyRepresentation {
    #[new]
    pub fn new(representation: String) -> Self {
        match representation.as_str() {
            "centroid" => PyRepresentation {
                inner: Representation::Centroid,
            },
            "profile" => PyRepresentation {
                inner: Representation::Profile,
            },
            _ => PyRepresentation {
                inner: Representation::Centroid,
            },
        }
    }

    #[getter]
    pub fn representation_as_string(&self) -> String {
        match self.inner {
            Representation::Centroid => "CENTROID",
            Representation::Profile => "PROFILE",
        }
        .to_string()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyProcessedSpectrum {
    pub inner: ProcessedSpectrum,
}

#[pymethods]
impl PyProcessedSpectrum {
    #[new]
    pub fn new(
        level: u8,
        id: String,
        file_id: usize,
        scan_start_time: f32,
        ion_injection_time: f32,
        precursors: Vec<PyPrecursor>,
        peaks: Vec<PyPeak>,
        total_ion_current: f32,
    ) -> Self {
        PyProcessedSpectrum {
            inner: ProcessedSpectrum {
                level,
                id,
                file_id,
                scan_start_time,
                ion_injection_time,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                peaks: peaks.into_iter().map(|p| p.inner).collect(),
                total_ion_current,
            },
        }
    }

    #[getter]
    pub fn level(&self) -> u8 {
        self.inner.level
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn scan_start_time(&self) -> f32 {
        self.inner.scan_start_time
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn precursors(&self) -> Vec<PyPrecursor> {
        self.inner
            .precursors
            .clone()
            .into_iter()
            .map(|p| PyPrecursor { inner: p })
            .collect()
    }

    #[getter]
    pub fn peaks(&self) -> Vec<PyPeak> {
        self.inner
            .peaks
            .clone()
            .into_iter()
            .map(|p| PyPeak { inner: p })
            .collect()
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }

    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        self.inner.extract_ms1_precursor()
    }

    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        self.inner.in_isolation_window(mz)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyRawSpectrum {
    pub inner: RawSpectrum,
}

#[pymethods]
impl PyRawSpectrum {
    #[new]
    pub fn new(
        file_id: usize,
        ms_level: u8,
        id: String,
        precursors: Vec<PyPrecursor>,
        representation: PyRepresentation,
        scan_start_time: f32,
        ion_injection_time: f32,
        total_ion_current: f32,
        mz: &PyArray1<f32>,
        intensity: &PyArray1<f32>,
    ) -> Self {
        let mz_vec = unsafe { mz.as_array().to_vec() };
        let intensity_vec = unsafe { intensity.as_array().to_vec() };

        PyRawSpectrum {
            inner: RawSpectrum {
                file_id,
                ms_level,
                id,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                representation: representation.inner,
                scan_start_time,
                ion_injection_time,
                total_ion_current,
                mz: mz_vec,
                intensity: intensity_vec,
            },
        }
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn ms_level(&self) -> u8 {
        self.inner.ms_level
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[getter]
    pub fn precursors(&self) -> Vec<PyPrecursor> {
        self.inner
            .precursors
            .clone()
            .into_iter()
            .map(|p| PyPrecursor { inner: p })
            .collect()
    }

    #[getter]
    pub fn representation(&self) -> PyRepresentation {
        PyRepresentation {
            inner: self.inner.representation,
        }
    }

    #[getter]
    pub fn scan_start_time(&self) -> f32 {
        self.inner.scan_start_time
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }

    #[getter]
    pub fn mz(&self, py: Python) -> Py<PyArray1<f32>> {
        self.inner.mz.clone().into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn intensity(&self, py: Python) -> Py<PyArray1<f32>> {
        self.inner.intensity.clone().into_pyarray(py).to_owned()
    }
}

#[pyclass]
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct PyPeak {
    pub inner: Peak,
}

#[pymethods]
impl PyPeak {
    #[new]
    pub fn new(mass: f32, intensity: f32) -> Self {
        PyPeak {
            inner: Peak { mass, intensity },
        }
    }

    #[getter]
    pub fn mass(&self) -> f32 {
        self.inner.mass
    }

    #[getter]
    pub fn intensity(&self) -> f32 {
        self.inner.intensity
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PySpectrumProcessor {
    pub inner: SpectrumProcessor,
}

#[pymethods]
impl PySpectrumProcessor {
    #[new]
    pub fn new(
        take_top_n: usize,
        max_fragment_mz: f32,
        min_fragment_mz: f32,
        deisotope: bool,
    ) -> Self {
        PySpectrumProcessor {
            inner: SpectrumProcessor {
                take_top_n,
                max_fragment_mz,
                min_fragment_mz,
                deisotope,
            },
        }
    }

    #[getter]
    pub fn take_top_n(&self) -> usize {
        self.inner.take_top_n
    }

    #[getter]
    pub fn max_fragment_mz(&self) -> f32 {
        self.inner.max_fragment_mz
    }

    #[getter]
    pub fn min_fragment_mz(&self) -> f32 {
        self.inner.min_fragment_mz
    }

    #[getter]
    pub fn deisotope(&self) -> bool {
        self.inner.deisotope
    }

    pub fn process(&self, spectrum: &PyRawSpectrum) -> PyProcessedSpectrum {
        PyProcessedSpectrum {
            inner: self.inner.process(spectrum.inner.clone()),
        }
    }
}

#[pyclass]
#[derive(PartialEq, Clone, Debug)]
pub struct PyDeisotoped {
    pub inner: Deisotoped,
}

#[pymethods]
impl PyDeisotoped {
    #[new]
    pub fn new(mz: f32, intensity: f32, charge: Option<u8>, envelope: Option<usize>) -> Self {
        PyDeisotoped {
            inner: Deisotoped {
                mz,
                intensity,
                charge,
                envelope,
            },
        }
    }

    #[getter]
    pub fn mz(&self) -> f32 {
        self.inner.mz
    }

    #[getter]
    pub fn intensity(&self) -> f32 {
        self.inner.intensity
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }

    #[getter]
    pub fn envelope(&self) -> Option<usize> {
        self.inner.envelope
    }
}

#[pyclass]
#[derive(Default, Clone, Debug)]
pub struct PyPrecursor {
    pub inner: Precursor,
}

#[pymethods]
impl PyPrecursor {
    #[new]
    pub fn new(
        mz: f32,
        intensity: Option<f32>,
        charge: Option<u8>,
        spectrum_ref: Option<String>,
        isolation_window: Option<PyTolerance>,
    ) -> Self {
        PyPrecursor {
            inner: Precursor {
                mz,
                intensity,
                charge,
                spectrum_ref,
                isolation_window: isolation_window.map(|t| t.inner),
            },
        }
    }

    #[getter]
    pub fn mz(&self) -> f32 {
        self.inner.mz
    }

    #[getter]
    pub fn intensity(&self) -> Option<f32> {
        self.inner.intensity
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }

    #[getter]
    pub fn spectrum_ref(&self) -> Option<String> {
        self.inner.spectrum_ref.clone()
    }

    #[getter]
    pub fn isolation_window(&self) -> Option<PyTolerance> {
        self.inner
            .isolation_window
            .clone()
            .map(|t| PyTolerance { inner: t })
    }
}

#[pymodule]
pub fn spectrum(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeak>()?;
    m.add_class::<PyDeisotoped>()?;
    m.add_class::<PyPrecursor>()?;
    m.add_class::<PySpectrumProcessor>()?;
    m.add_class::<PyRepresentation>()?;
    m.add_class::<PyRawSpectrum>()?;
    m.add_class::<PyProcessedSpectrum>()?;
    Ok(())
}
