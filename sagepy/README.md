# SAGEpy
A python interface to the core [SAGE](https://github.com/lazear/sage) search engine for mass spectrometry proteomics

<p align="center">
  <img src="sagepy_logo.png" alt="logo" width="250"/>
</p>

## Installation
`sagepy` is now available via pip:
```
pip install sagepy
```


### Build from source

1. Clone our fork of the SAGE repository:
```
git clone git@github.com:theGreatHerrLebert/sage.git
```

2. Install the sage-core bindings using maturin, optionally in a virtual environment:
```
cd sage/crates/sagepy-connector

# Install maturin
pip install maturin

# Build and install the bindings
maturin build --release

# Install the bindings
pip install target/wheels/sagepy_connector-0.1.0-cp38-cp38-manylinux2014_x86_64.whl [--force-reinstall]
```
This will provide you with a python exposed version of the core SAGE library.

3. Install the sagepy python package with poetry:
```
git clone git@github.com:theGreatHerrLebert/sagepy.git

cd sagepy

# Install poetry
pip install poetry

# Install sagepy
poetry install
```

## Usage
`sagepy` is a python interface to the core SAGE search engine. It exposes
the core functionality of SAGE in a pythonic way, allowing you to use it for a direct integration 
into your python-based proteomics workflow. So far, it mainly mirrors structs that are available
in the core SAGE library. 

### Example generation of a sage database
```python
import numpy as np
from sagepy.core import EnzymeBuilder, SageSearchConfiguration

# configure a trypsin-like digestor of fasta files
enzyme_builder = EnzymeBuilder(
    missed_cleavages=2, 
    min_len=5, 
    max_len=50, 
    cleave_at='KR', 
    restrict='P', 
    c_terminal=True,
)

# UPDATE: Modification handling is simplified, using canonical UNIMOD notation
static_mods = {"C": "[UNIMOD:4]"}  # static cysteine modification
variable_mods = {"M": ["[UNIMOD:35]"]}

with open('path/to/reference.fasta', 'r') as infile:
    fasta = infile.read()

# set-up a config for a sage-database
sage_config = SageSearchConfiguration(
    fasta=fasta,
    static_mods=static_mods,
    variable_mods=variable_mods,
    enzyme_builder=enzyme_builder,
    generate_decoys=True,
    bucket_size=int(np.power(2, 14))
)

# generate the database for searching against
indexed_db = sage_config.generate_indexed_database()
```

### Generate a query
```python
import numpy as np
from sagepy.core import Precursor, RawSpectrum, ProcessedSpectrum, SpectrumProcessor, Tolerance, Scorer, Representation

### Example search of a sage database
precursor = Precursor(
    charge=2,
    mz=506.77,
)

intensity = np.array([ 202.,  170.,  205.,  152., 1069.,  595.,  198.,  805.,  187.,
        194.,  197.,  169.,  196.,  209.,  638.,  372.,  235.,  399.,
        194.,  185.,  181.,  170.,  407.,  150.,  157.,  175.,  273.,
       1135.,  881.,  337.,  311.,  243.,  310.,  153.,  162.,  210.,
        277.,  206.,  189.,  259.,  658.,  383.,  166.,  169.,  219.,
        186.,  221.,  193.,  367.,  283.,  237.,  157.,  372., 1276.,
       1618., 1102.,  404.,  232.,  456.,  765.,  507.,  223.,  258.,
        402.,  187.,  158.,  153.,  304.,  218.,  223.,  156., 1605.,
       1165., 1062.,  434.,  208.,  155.,  197.,  221.,  697.,  397.,
        180.,  195.,  512.,  252.,  367.,  305.,  335.,  175.,  174.,
        296.,  212.], dtype=np.float32)

mz = np.array([272.16873692, 356.16844797, 406.71079396, 406.71396814,
       406.71714233, 406.72031653, 407.21246768, 407.21564382,
       407.21881996, 407.22199612, 407.7144506 , 407.71762869,
       488.27537883, 488.28581266, 499.29228981, 499.29580676,
       499.29932372, 499.30284069, 506.75478369, 507.26157767,
       541.26272227, 553.29188809, 577.30432041, 577.30810217,
       595.32672633, 597.2907525 , 603.27568881, 614.32036769,
       614.32426881, 614.32816995, 615.3272682 , 615.33117252,
       616.33108578, 617.33572156, 636.30924838, 637.30619081,
       637.31016425, 665.36284673, 666.36197292, 674.35335834,
       674.35744565, 674.36153297, 675.35511968, 675.36330039,
       679.3531909 , 680.35044702, 680.35455247, 687.36822726,
       687.37648041, 688.37547678, 697.3616813 , 700.3617026 ,
       715.36157366, 715.36578342, 715.36999319, 715.37420297,
       715.37841277, 715.38262258, 716.36384605, 716.37227148,
       716.38069696, 717.37103577, 725.35228543, 749.39291293,
       749.39722166, 750.38424802, 786.44692356, 786.45575152,
       787.4492132 , 787.45804678, 795.39284711, 812.41777208,
       812.42225834, 812.42674462, 812.4312309 , 812.44020351,
       813.40504794, 813.41851494, 813.42300396, 813.427493  ,
       813.43198205, 813.44544927, 814.43784098, 828.42202737,
       828.4265576 , 851.43464868, 899.45327427, 899.46271517,
       912.45278821, 913.44673363, 915.45053417, 915.46482091], dtype=np.float32)

raw_spectrum = RawSpectrum(
    file_id=1,
    spec_id='DEMO-SPEC',
    total_ion_current=12667.0,
    precursors=[precursor],
    mz=mz,
    intensity=intensity
)

spec_processor = SpectrumProcessor(take_top_n=75)
query = spec_processor.process(raw_spectrum)
```

### Search a database
```python
from sagepy.core import Scorer

# UPDATE: pass modifications to the scorer, necessary for PTM handling
scorer = Scorer(report_psms=2, min_matched_peaks=5, variable_mods=variable_mods, static_mods=static_mods)
results = scorer.score(db=indexed_db, spectrum=query)
```

potential output:
```
[Feature(idx: PeptideIx(1009105), peptide_len: 9, spec_id: DEMO-SPEC, file_id: 1, rank: 1, label: 1, exp. mass: 1011.5254516601562, cal. mass: 1011.5347900390625, charge: 2, retention time: 0.0, aligned rt: 0.0, predicted rt: 0.0, delta rt model: 0.9990000128746033, delta mass: 2989.41943359375, isotope error: 3.010050058364868, average ppm: 5.889466285705566, hyperscore: 15.020833459653923, delta_next: 0.0, delta_best: 0.0, matched peaks: 5, longest b: 0,longest y: 4, longest y pct: 0.4444444477558136, missed cleavages: 0, matched intensity pct: 14.81151294708252, scored candidates: 9340, poisson: -2.177020383746938, discriminant score: 0.0, posterior error: 1.0, spectrum q: 1.0, peptide q: 1.0, protein q: 1.0, ms2 intensity: 4652.0, ms1 intensity: 0.0), Feature(idx: PeptideIx(1009105), peptide_len: 9, spec_id: DEMO-SPEC, file_id: 1, rank: 2, label: 1, exp. mass: 1011.5254516601562, cal. mass: 1011.5347900390625, charge: 2, retention time: 0.0, aligned rt: 0.0, predicted rt: 0.0, delta rt model: 0.9990000128746033, delta mass: 1001.641845703125, isotope error: 1.003350019454956, average ppm: 5.889466285705566, hyperscore: 15.020833459653923, delta_next: 0.0, delta_best: 0.0, matched peaks: 5, longest b: 0,longest y: 4, longest y pct: 0.4444444477558136, missed cleavages: 0, matched intensity pct: 14.81151294708252, scored candidates: 9340, poisson: -2.177020383746938, discriminant score: 0.0, posterior error: 1.0, spectrum q: 1.0, peptide q: 1.0, protein q: 1.0, ms2 intensity: 4652.0, ms1 intensity: 0.0)]
```
