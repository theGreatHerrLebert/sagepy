use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::prelude::*;

use crate::py_mass::PyTolerance;
use sage_core::spectrum::{Deisotoped, IMPeak, Peak, Precursor, ProcessedSpectrum, RawSpectrum, Representation, SpectrumProcessor};

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
    pub inner: ProcessedSpectrum<Peak>,
    pub collision_energies: Vec<Option<f32>>,
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
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();
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
            collision_energies,
        }
    }

    #[getter]
    pub fn level(&self) -> u8 {
        self.inner.level
    }

    #[setter]
    pub fn set_level(&mut self, level: u8) {
        self.inner.level = level;
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[setter]
    pub fn set_id(&mut self, id: String) {
        self.inner.id = id;
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[setter]
    pub fn set_file_id(&mut self, file_id: usize) {
        self.inner.file_id = file_id;
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
            .precursors.iter().zip(self.collision_energies.iter())
            .into_iter()
            .map(|p| PyPrecursor { inner: p.0.clone(), collision_energy: *p.1})
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

    pub fn calibrate_mz_ppm(&mut self, ppm: f32) {
        for peak in self.inner.peaks.iter_mut() {
            let ppm_mass = peak.mass / 1e6;
            let ppm_error = ppm_mass * ppm;
            peak.mass -= ppm_error;
        }
    }

    #[getter]
    pub fn collision_energies(&self) -> Vec<Option<f32>> {
        self.collision_energies.clone()
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
pub struct PyProcessedIMSpectrum {
    pub inner: ProcessedSpectrum<IMPeak>,
    pub collision_energies: Vec<Option<f32>>,
}

#[pymethods]
impl PyProcessedIMSpectrum {
    #[new]
    pub fn new(
        level: u8,
        id: String,
        file_id: usize,
        scan_start_time: f32,
        ion_injection_time: f32,
        precursors: Vec<PyPrecursor>,
        peaks: Vec<PyIMPeak>,
        total_ion_current: f32,
    ) -> Self {
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();

        PyProcessedIMSpectrum {
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
            collision_energies,
        }
    }

    #[getter]
    pub fn level(&self) -> u8 {
        self.inner.level
    }

    #[setter]
    pub fn set_level(&mut self, level: u8) {
        self.inner.level = level;
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[setter]
    pub fn set_id(&mut self, id: String) {
        self.inner.id = id;
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[setter]
    pub fn set_file_id(&mut self, file_id: usize) {
        self.inner.file_id = file_id;
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
            .iter()
            .zip(self.collision_energies.iter())
            .map(|(p, ce)| PyPrecursor {
                inner: p.clone(),
                collision_energy: *ce,
            })
            .collect()
    }

    #[getter]
    pub fn peaks(&self) -> Vec<PyIMPeak> {
        self.inner
            .peaks
            .clone()
            .into_iter()
            .map(|p| PyIMPeak { inner: p })
            .collect()
    }

    #[getter]
    pub fn collision_energies(&self) -> Vec<Option<f32>> {
        self.collision_energies.clone()
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyRawSpectrum {
    pub inner: RawSpectrum,
    pub collision_energies: Vec<Option<f32>>,
}

#[pymethods]
impl PyRawSpectrum {
    #[new]
    #[pyo3(signature = (file_id, ms_level, id, precursors, representation, scan_start_time, ion_injection_time, total_ion_current, mz, intensity, mobility=None))]
    pub fn new(
        file_id: usize,
        ms_level: u8,
        id: String,
        precursors: Vec<PyPrecursor>,
        representation: PyRepresentation,
        scan_start_time: f32,
        ion_injection_time: f32,
        total_ion_current: f32,
        mz: &Bound<'_, PyArray1<f32>>,
        intensity: &Bound<'_, PyArray1<f32>>,
        mobility: Option<Vec<f32>>,
    ) -> Self {
        let mz_vec = unsafe { mz.as_array().to_vec() };
        let intensity_vec = unsafe { intensity.as_array().to_vec() };
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();

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
                mobility,
            },
            collision_energies,
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
            .precursors.iter().zip(self.collision_energies.iter())
            .into_iter()
            .map(|p| PyPrecursor { inner: p.0.clone(), collision_energy: *p.1 })
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
        self.inner.mz.clone().into_pyarray(py).unbind()
    }

    #[getter]
    pub fn intensity(&self, py: Python) -> Py<PyArray1<f32>> {
        self.inner.intensity.clone().into_pyarray(py).unbind()
    }
    
    #[getter]
    pub fn mobility(&self, py: Python) -> Option<Py<PyArray1<f32>>> {
        self.inner.mobility.as_ref().map(|mob| mob.clone().into_pyarray(py).unbind())
    }

    pub fn filter_top_n(&self, n: usize) -> PyRawSpectrum {
        if let Some(mobility_vec) = &self.inner.mobility {
            // Define inline struct for 3D peak
            struct Peak3D {
                mass: f32,
                intensity: f32,
                mobility: f32,
            }

            let mz_iter = self.inner.mz.iter();
            let intensity_iter = self.inner.intensity.iter();

            let mut peaks = mz_iter
                .zip(intensity_iter)
                .zip(mobility_vec.iter())
                .map(|((m, i), mob)| Peak3D {
                    mass: *m,
                    intensity: *i,
                    mobility: *mob,
                })
                .collect::<Vec<_>>();

            peaks.sort_by(|a, b| b.intensity.partial_cmp(&a.intensity).unwrap_or(std::cmp::Ordering::Equal));
            let mut top_peaks = peaks.into_iter().take(n).collect::<Vec<_>>();
            top_peaks.sort_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap_or(std::cmp::Ordering::Equal));

            PyRawSpectrum {
                inner: RawSpectrum {
                    file_id: self.inner.file_id,
                    ms_level: self.inner.ms_level,
                    id: self.inner.id.clone(),
                    precursors: self.inner.precursors.clone(),
                    representation: self.inner.representation,
                    scan_start_time: self.inner.scan_start_time,
                    ion_injection_time: self.inner.ion_injection_time,
                    total_ion_current: self.inner.total_ion_current,
                    mz: top_peaks.iter().map(|p| p.mass).collect(),
                    intensity: top_peaks.iter().map(|p| p.intensity).collect(),
                    mobility: Some(top_peaks.iter().map(|p| p.mobility).collect()),
                },
                collision_energies: self.collision_energies.clone(),
            }
        } else {
            // Define inline struct for 2D peak
            struct Peak {
                mass: f32,
                intensity: f32,
            }

            let mut peaks = self.inner.mz.iter()
                .zip(self.inner.intensity.iter())
                .map(|(m, i)| Peak { mass: *m, intensity: *i })
                .collect::<Vec<_>>();

            peaks.sort_by(|a, b| b.intensity.partial_cmp(&a.intensity).unwrap_or(std::cmp::Ordering::Equal));
            let mut top_peaks = peaks.into_iter().take(n).collect::<Vec<_>>();
            top_peaks.sort_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap_or(std::cmp::Ordering::Equal));

            PyRawSpectrum {
                inner: RawSpectrum {
                    file_id: self.inner.file_id,
                    ms_level: self.inner.ms_level,
                    id: self.inner.id.clone(),
                    precursors: self.inner.precursors.clone(),
                    representation: self.inner.representation,
                    scan_start_time: self.inner.scan_start_time,
                    ion_injection_time: self.inner.ion_injection_time,
                    total_ion_current: self.inner.total_ion_current,
                    mz: top_peaks.iter().map(|p| p.mass).collect(),
                    intensity: top_peaks.iter().map(|p| p.intensity).collect(),
                    mobility: None,
                },
                collision_energies: self.collision_energies.clone(),
            }
        }
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
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct PyIMPeak {
    pub inner: IMPeak,
}

#[pymethods]
impl PyIMPeak {
    #[new]
    pub fn new(mass: f32, intensity: f32, mobility: f32) -> Self {
        PyIMPeak {
            inner: IMPeak {
                mass,
                intensity,
                mobility,
            },
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
    
    #[getter]
    pub fn mobility(&self) -> f32 {
        self.inner.mobility
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
        min_deisotope_mz: f32,
        deisotope: bool,
    ) -> Self {
        PySpectrumProcessor {
            inner: SpectrumProcessor {
                take_top_n,
                min_deisotope_mz,
                deisotope,
            },
        }
    }

    #[getter]
    pub fn take_top_n(&self) -> usize {
        self.inner.take_top_n
    }

    #[getter]
    pub fn deisotope(&self) -> bool {
        self.inner.deisotope
    }

    pub fn process(&self, spectrum: &PyRawSpectrum) -> PyProcessedSpectrum {
        PyProcessedSpectrum {
            inner: self.inner.process(spectrum.inner.clone()),
            collision_energies: spectrum.collision_energies.clone(),
        }
    }
    pub fn process_with_mobility(&self, spectrum: &PyRawSpectrum) -> PyProcessedIMSpectrum {
        PyProcessedIMSpectrum {
            inner: self.inner.process_with_mobility(spectrum.inner.clone()),
            collision_energies: spectrum.collision_energies.clone(),
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
    #[pyo3(signature = (mz, intensity, charge=None, envelope=None))]
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
    pub collision_energy: Option<f32>,
}

#[pymethods]
impl PyPrecursor {
    #[new]
    #[pyo3(signature = (mz, intensity=None, charge=None, spectrum_ref=None, isolation_window=None, inverse_ion_mobility=None, collision_energy=None))]
    pub fn new(
        mz: f32,
        intensity: Option<f32>,
        charge: Option<u8>,
        spectrum_ref: Option<String>,
        isolation_window: Option<PyTolerance>,
        inverse_ion_mobility: Option<f32>,
        collision_energy: Option<f32>,
    ) -> Self {
        PyPrecursor {
            inner: Precursor {
                mz,
                intensity,
                charge,
                spectrum_ref,
                isolation_window: isolation_window.map(|t| t.inner),
                inverse_ion_mobility,
            },
            collision_energy,
        }
    }

    #[getter]
    pub fn mz(&self) -> f32 {
        self.inner.mz
    }

    pub fn calibrate_mz_ppm(&mut self, ppm: f32) {
        let mz_ppm = self.inner.mz / 1e6;
        let ppm_error = mz_ppm * ppm;
        self.inner.mz -= ppm_error;
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

    #[getter]
    pub fn inverse_ion_mobility(&self) -> Option<f32> {
        self.inner.inverse_ion_mobility
    }

    #[getter]
    pub fn collision_energy(&self) -> Option<f32> {
        self.collision_energy
    }
}

#[pymodule]
pub fn py_spectrum(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPeak>()?;
    m.add_class::<PyIMPeak>()?;
    m.add_class::<PyDeisotoped>()?;
    m.add_class::<PyPrecursor>()?;
    m.add_class::<PySpectrumProcessor>()?;
    m.add_class::<PyRepresentation>()?;
    m.add_class::<PyRawSpectrum>()?;
    m.add_class::<PyProcessedSpectrum>()?;
    m.add_class::<PyProcessedIMSpectrum>()?;
    Ok(())
}
