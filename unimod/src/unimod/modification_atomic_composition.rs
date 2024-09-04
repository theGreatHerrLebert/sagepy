use std::collections::HashMap;

/// Unimod Modifications
///
/// # Arguments
///
/// None
///
/// # Returns
///
/// * `HashMap<String, HashMap<&'static str, i32>>` - a map of unimod modification names to their atomic compositions
///
/// # Example
///
/// ```
/// use unimod::unimod::modification_atomic_composition;
/// use std::collections::HashMap;
///
/// let composition = modification_atomic_composition();
/// assert_eq!(composition.get("[UNIMOD:1]"), Some(&HashMap::from([("C", 2), ("H", 2), ("O", 1)])));
/// ```
pub fn modification_atomic_composition() -> HashMap<String, HashMap<&'static str, i32>> {
    let mut composition: HashMap<String, HashMap<&'static str, i32>> = HashMap::new();
    composition.insert("[UNIMOD:1]".to_string(), HashMap::from([("H", 2), ("C", 2), ("O", 1)])); // Acetyl
    composition.insert("[UNIMOD:2]".to_string(), HashMap::from([("H", 1), ("N", 1), ("O", -1)])); // Amidated
    composition.insert("[UNIMOD:3]".to_string(), HashMap::from([("H", 14), ("C", 10), ("N", 2), ("O", 2), ("S", 1)])); // Biotin
    composition.insert("[UNIMOD:4]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1)])); // Carbamidomethyl
    composition.insert("[UNIMOD:5]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", 1), ("O", 1)])); // Carbamyl
    composition.insert("[UNIMOD:6]".to_string(), HashMap::from([("H", 2), ("C", 2), ("O", 2)])); // Carboxymethyl
    composition.insert("[UNIMOD:7]".to_string(), HashMap::from([("H", -1), ("N", -1), ("O", 1)])); // Deamidated
    composition.insert("[UNIMOD:8]".to_string(), HashMap::from([("H", 38), ("C", 22), ("N", 4), ("O", 6), ("S", 1)])); // ICAT-G
    composition.insert("[UNIMOD:9]".to_string(), HashMap::from([("H", 30), ("2H", 8), ("C", 22), ("N", 4), ("O", 6), ("S", 1)])); // ICAT-G:2H(8)
    composition.insert("[UNIMOD:10]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", 1), ("S", -1)])); // Met->Hse
    composition.insert("[UNIMOD:11]".to_string(), HashMap::from([("H", -4), ("C", -1), ("S", -1)])); // Met->Hsl
    composition.insert("[UNIMOD:12]".to_string(), HashMap::from([("H", 26), ("2H", 8), ("C", 20), ("N", 4), ("O", 5), ("S", 1)])); // ICAT-D:2H(8)
    composition.insert("[UNIMOD:13]".to_string(), HashMap::from([("H", 34), ("C", 20), ("N", 4), ("O", 5), ("S", 1)])); // ICAT-D
    composition.insert("[UNIMOD:17]".to_string(), HashMap::from([("H", 9), ("C", 5), ("N", 1), ("O", 1)])); // NIPCAM
    composition.insert("[UNIMOD:20]".to_string(), HashMap::from([("H", 30), ("C", 18), ("N", 4), ("O", 5), ("S", 1)])); // PEO-Iodoacetyl-LC-Biotin
    composition.insert("[UNIMOD:21]".to_string(), HashMap::from([("H", 1), ("O", 3), ("P", 1)])); // Phospho
    composition.insert("[UNIMOD:23]".to_string(), HashMap::from([("H", -2), ("O", -1)])); // Dehydrated
    composition.insert("[UNIMOD:24]".to_string(), HashMap::from([("H", 5), ("C", 3), ("N", 1), ("O", 1)])); // Propionamide
    composition.insert("[UNIMOD:25]".to_string(), HashMap::from([("H", 5), ("C", 7), ("N", 1), ("O", 1)])); // Pyridylacetyl
    composition.insert("[UNIMOD:26]".to_string(), HashMap::from([("C", 2), ("O", 1)])); // Pyro-carbamidomethyl
    composition.insert("[UNIMOD:27]".to_string(), HashMap::from([("H", -2), ("O", -1)])); // Glu->pyro-Glu
    composition.insert("[UNIMOD:28]".to_string(), HashMap::from([("H", -3), ("N", -1)])); // Gln->pyro-Glu
    composition.insert("[UNIMOD:29]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 2)])); // SMA
    composition.insert("[UNIMOD:30]".to_string(), HashMap::from([("H", -1), ("Na", 1)])); // Cation:Na
    composition.insert("[UNIMOD:31]".to_string(), HashMap::from([("H", 7), ("C", 7), ("N", 1)])); // Pyridylethyl
    composition.insert("[UNIMOD:34]".to_string(), HashMap::from([("H", 2), ("C", 1)])); // Methyl
    composition.insert("[UNIMOD:35]".to_string(), HashMap::from([("O", 1)])); // Oxidation
    composition.insert("[UNIMOD:36]".to_string(), HashMap::from([("H", 4), ("C", 2)])); // Dimethyl
    composition.insert("[UNIMOD:37]".to_string(), HashMap::from([("H", 6), ("C", 3)])); // Trimethyl
    composition.insert("[UNIMOD:39]".to_string(), HashMap::from([("H", 2), ("C", 1), ("S", 1)])); // Methylthio
    composition.insert("[UNIMOD:40]".to_string(), HashMap::from([("O", 3), ("S", 1)])); // Sulfo
    composition.insert("[UNIMOD:41]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 5)])); // Hex
    composition.insert("[UNIMOD:42]".to_string(), HashMap::from([("H", 12), ("C", 8), ("O", 1), ("S", 2)])); // Lipoyl
    composition.insert("[UNIMOD:43]".to_string(), HashMap::from([("H", 13), ("C", 8), ("N", 1), ("O", 5)])); // HexNAc
    composition.insert("[UNIMOD:44]".to_string(), HashMap::from([("H", 24), ("C", 15)])); // Farnesyl
    composition.insert("[UNIMOD:45]".to_string(), HashMap::from([("H", 26), ("C", 14), ("O", 1)])); // Myristoyl
    composition.insert("[UNIMOD:46]".to_string(), HashMap::from([("H", 8), ("C", 8), ("N", 1), ("O", 5), ("P", 1)])); // PyridoxalPhosphate
    composition.insert("[UNIMOD:47]".to_string(), HashMap::from([("H", 30), ("C", 16), ("O", 1)])); // Palmitoyl
    composition.insert("[UNIMOD:48]".to_string(), HashMap::from([("H", 32), ("C", 20)])); // GeranylGeranyl
    composition.insert("[UNIMOD:49]".to_string(), HashMap::from([("H", 21), ("C", 11), ("N", 2), ("O", 6), ("P", 1), ("S", 1)])); // Phosphopantetheine
    composition.insert("[UNIMOD:50]".to_string(), HashMap::from([("H", 31), ("C", 27), ("N", 9), ("O", 15), ("P", 2)])); // FAD
    composition.insert("[UNIMOD:51]".to_string(), HashMap::from([("H", 96), ("C", 51), ("O", 5)])); // Tripalmitate
    composition.insert("[UNIMOD:52]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 2)])); // Guanidinyl
    composition.insert("[UNIMOD:53]".to_string(), HashMap::from([("H", 16), ("C", 9), ("O", 2)])); // HNE
    composition.insert("[UNIMOD:54]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 6)])); // Glucuronyl
    composition.insert("[UNIMOD:55]".to_string(), HashMap::from([("H", 15), ("C", 10), ("N", 3), ("O", 6), ("S", 1)])); // Glutathione
    composition.insert("[UNIMOD:56]".to_string(), HashMap::from([("H", -1), ("2H", 3), ("C", 2), ("O", 1)])); // Acetyl:2H(3)
    composition.insert("[UNIMOD:58]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 1)])); // Propionyl
    composition.insert("[UNIMOD:59]".to_string(), HashMap::from([("H", 4), ("13C", 3), ("O", 1)])); // Propionyl:13C(3)
    composition.insert("[UNIMOD:60]".to_string(), HashMap::from([("H", 13), ("C", 7), ("N", 1), ("O", 1)])); // GIST-Quat
    composition.insert("[UNIMOD:61]".to_string(), HashMap::from([("H", 10), ("2H", 3), ("C", 7), ("N", 1), ("O", 1)])); // GIST-Quat:2H(3)
    composition.insert("[UNIMOD:62]".to_string(), HashMap::from([("H", 7), ("2H", 6), ("C", 7), ("N", 1), ("O", 1)])); // GIST-Quat:2H(6)
    composition.insert("[UNIMOD:63]".to_string(), HashMap::from([("H", 4), ("2H", 9), ("C", 7), ("N", 1), ("O", 1)])); // GIST-Quat:2H(9)
    composition.insert("[UNIMOD:64]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 3)])); // Succinyl
    composition.insert("[UNIMOD:65]".to_string(), HashMap::from([("2H", 4), ("C", 4), ("O", 3)])); // Succinyl:2H(4)
    composition.insert("[UNIMOD:66]".to_string(), HashMap::from([("H", 4), ("13C", 4), ("O", 3)])); // Succinyl:13C(4)
    composition.insert("[UNIMOD:89]".to_string(), HashMap::from([("H", 15), ("C", 10), ("N", 3), ("O", 1), ("S", 1)])); // Iminobiotin
    composition.insert("[UNIMOD:90]".to_string(), HashMap::from([("H", 26), ("C", 16), ("N", 4), ("O", 2), ("S", 1)])); // ESP
    composition.insert("[UNIMOD:91]".to_string(), HashMap::from([("H", 16), ("2H", 10), ("C", 16), ("N", 4), ("O", 2), ("S", 1)])); // ESP:2H(10)
    composition.insert("[UNIMOD:92]".to_string(), HashMap::from([("H", 25), ("C", 16), ("N", 3), ("O", 3), ("S", 1)])); // NHS-LC-Biotin
    composition.insert("[UNIMOD:93]".to_string(), HashMap::from([("H", 39), ("C", 25), ("N", 5), ("O", 6), ("S", 3)])); // EDT-maleimide-PEO-biotin
    composition.insert("[UNIMOD:94]".to_string(), HashMap::from([("H", 4), ("C", 3), ("N", 2)])); // IMID
    composition.insert("[UNIMOD:95]".to_string(), HashMap::from([("2H", 4), ("C", 3), ("N", 2)])); // IMID:2H(4)
    composition.insert("[UNIMOD:97]".to_string(), HashMap::from([("H", 2), ("2H", 3), ("C", 3), ("N", 1), ("O", 1)])); // Propionamide:2H(3)
    composition.insert("[UNIMOD:105]".to_string(), HashMap::from([("H", 17), ("C", 10), ("N", 3), ("O", 3)])); // ICAT-C
    composition.insert("[UNIMOD:106]".to_string(), HashMap::from([("H", 17), ("C", 1), ("13C", 9), ("N", 3), ("O", 3)])); // ICAT-C:13C(9)
    composition.insert("[UNIMOD:107]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 2), ("S", 1)])); // FormylMet
    composition.insert("[UNIMOD:108]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 1), ("O", 2)])); // Nethylmaleimide
    composition.insert("[UNIMOD:112]".to_string(), HashMap::from([("H", 26), ("C", 16), ("N", 4), ("O", 3), ("S", 1)])); // OxLysBiotinRed
    composition.insert("[UNIMOD:113]".to_string(), HashMap::from([("H", 24), ("C", 16), ("N", 4), ("O", 3), ("S", 1)])); // OxLysBiotin
    composition.insert("[UNIMOD:114]".to_string(), HashMap::from([("H", 29), ("C", 16), ("N", 5), ("O", 3), ("S", 1)])); // OxProBiotinRed
    composition.insert("[UNIMOD:115]".to_string(), HashMap::from([("H", 27), ("C", 16), ("N", 5), ("O", 3), ("S", 1)])); // OxProBiotin
    composition.insert("[UNIMOD:116]".to_string(), HashMap::from([("H", 22), ("C", 15), ("N", 2), ("O", 3), ("S", 1)])); // OxArgBiotin
    composition.insert("[UNIMOD:117]".to_string(), HashMap::from([("H", 24), ("C", 15), ("N", 2), ("O", 3), ("S", 1)])); // OxArgBiotinRed
    composition.insert("[UNIMOD:118]".to_string(), HashMap::from([("H", 34), ("C", 20), ("N", 4), ("O", 4), ("S", 3)])); // EDT-iodoacetyl-PEO-biotin
    composition.insert("[UNIMOD:119]".to_string(), HashMap::from([("H", 21), ("C", 22), ("P", 1)])); // IBTP
    composition.insert("[UNIMOD:121]".to_string(), HashMap::from([("H", 6), ("C", 4), ("N", 2), ("O", 2)])); // GG
    composition.insert("[UNIMOD:122]".to_string(), HashMap::from([("C", 1), ("O", 1)])); // Formyl
    composition.insert("[UNIMOD:123]".to_string(), HashMap::from([("H", 20), ("C", 15), ("N", 1), ("O", 6), ("Cl", 1)])); // ICAT-H
    composition.insert("[UNIMOD:124]".to_string(), HashMap::from([("H", 20), ("C", 9), ("13C", 6), ("N", 1), ("O", 6), ("Cl", 1)])); // ICAT-H:13C(6)
    composition.insert("[UNIMOD:126]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 1), ("S", 1)])); // Xlink:DTSSP[88]
    composition.insert("[UNIMOD:127]".to_string(), HashMap::from([("H", -1), ("F", 1)])); // Fluoro
    composition.insert("[UNIMOD:128]".to_string(), HashMap::from([("H", 13), ("C", 22), ("N", 1), ("O", 6)])); // Fluorescein
    composition.insert("[UNIMOD:129]".to_string(), HashMap::from([("H", -1), ("I", 1)])); // Iodo
    composition.insert("[UNIMOD:130]".to_string(), HashMap::from([("H", -2), ("I", 2)])); // Diiodo
    composition.insert("[UNIMOD:131]".to_string(), HashMap::from([("H", -3), ("I", 3)])); // Triiodo
    composition.insert("[UNIMOD:134]".to_string(), HashMap::from([("H", 24), ("C", 14), ("O", 1)])); // Myristoleyl
    composition.insert("[UNIMOD:135]".to_string(), HashMap::from([("H", 22), ("C", 14), ("O", 1)])); // Myristoyl+Delta:H(-4)
    composition.insert("[UNIMOD:136]".to_string(), HashMap::from([("H", 4), ("C", 7), ("O", 1)])); // Benzoyl
    composition.insert("[UNIMOD:137]".to_string(), HashMap::from([("H", 76), ("C", 46), ("N", 2), ("O", 35)])); // Hex(5)HexNAc(2)
    composition.insert("[UNIMOD:139]".to_string(), HashMap::from([("H", 11), ("C", 12), ("N", 1), ("O", 2), ("S", 1)])); // Dansyl
    composition.insert("[UNIMOD:140]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", -2)])); // a-type-ion
    composition.insert("[UNIMOD:141]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1)])); // Amidine
    composition.insert("[UNIMOD:142]".to_string(), HashMap::from([("H", 23), ("C", 14), ("N", 1), ("O", 9)])); // HexNAc(1)dHex(1)
    composition.insert("[UNIMOD:143]".to_string(), HashMap::from([("H", 26), ("C", 16), ("N", 2), ("O", 10)])); // HexNAc(2)
    composition.insert("[UNIMOD:144]".to_string(), HashMap::from([("H", 30), ("C", 18), ("O", 15)])); // Hex(3)
    composition.insert("[UNIMOD:145]".to_string(), HashMap::from([("H", 33), ("C", 20), ("N", 1), ("O", 13)])); // HexNAc(1)dHex(2)
    composition.insert("[UNIMOD:146]".to_string(), HashMap::from([("H", 33), ("C", 20), ("N", 1), ("O", 14)])); // Hex(1)HexNAc(1)dHex(1)
    composition.insert("[UNIMOD:147]".to_string(), HashMap::from([("H", 36), ("C", 22), ("N", 2), ("O", 14)])); // HexNAc(2)dHex(1)
    composition.insert("[UNIMOD:148]".to_string(), HashMap::from([("H", 36), ("C", 22), ("N", 2), ("O", 15)])); // Hex(1)HexNAc(2)
    composition.insert("[UNIMOD:149]".to_string(), HashMap::from([("H", 40), ("C", 25), ("N", 2), ("O", 18)])); // Hex(1)HexNAc(1)NeuAc(1)
    composition.insert("[UNIMOD:150]".to_string(), HashMap::from([("H", 46), ("C", 28), ("N", 2), ("O", 18)])); // HexNAc(2)dHex(2)
    composition.insert("[UNIMOD:151]".to_string(), HashMap::from([("H", 44), ("C", 27), ("N", 2), ("O", 19)])); // Hex(1)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:152]".to_string(), HashMap::from([("H", 46), ("C", 28), ("N", 2), ("O", 19)])); // Hex(1)HexNAc(2)dHex(1)
    composition.insert("[UNIMOD:153]".to_string(), HashMap::from([("H", 46), ("C", 28), ("N", 2), ("O", 20)])); // Hex(2)HexNAc(2)
    composition.insert("[UNIMOD:154]".to_string(), HashMap::from([("H", 51), ("C", 31), ("N", 1), ("O", 24)])); // Hex(3)HexNAc(1)Pent(1)
    composition.insert("[UNIMOD:155]".to_string(), HashMap::from([("H", 54), ("C", 33), ("N", 2), ("O", 23)])); // Hex(1)HexNAc(2)dHex(1)Pent(1)
    composition.insert("[UNIMOD:156]".to_string(), HashMap::from([("H", 56), ("C", 34), ("N", 2), ("O", 23)])); // Hex(1)HexNAc(2)dHex(2)
    composition.insert("[UNIMOD:157]".to_string(), HashMap::from([("H", 54), ("C", 33), ("N", 2), ("O", 24)])); // Hex(2)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:158]".to_string(), HashMap::from([("H", 56), ("C", 34), ("N", 2), ("O", 24)])); // Hex(2)HexNAc(2)dHex(1)
    composition.insert("[UNIMOD:159]".to_string(), HashMap::from([("H", 56), ("C", 34), ("N", 2), ("O", 25)])); // Hex(3)HexNAc(2)
    composition.insert("[UNIMOD:160]".to_string(), HashMap::from([("H", 57), ("C", 36), ("N", 3), ("O", 26)])); // Hex(1)HexNAc(1)NeuAc(2)
    composition.insert("[UNIMOD:161]".to_string(), HashMap::from([("H", 57), ("C", 34), ("N", 2), ("O", 28), ("P", 1)])); // Hex(3)HexNAc(2)Phos(1)
    composition.insert("[UNIMOD:162]".to_string(), HashMap::from([("S", -1), ("Se", 1)])); // Delta:S(-1)Se(1)
    composition.insert("[UNIMOD:170]".to_string(), HashMap::from([("H", -1), ("N", -1), ("18O", 1)])); // Delta:H(-1)N(-1)18O(1)
    composition.insert("[UNIMOD:171]".to_string(), HashMap::from([("H", 3), ("13C", 6), ("N", 1), ("O", 2), ("S", 1)])); // NBS:13C(6)
    composition.insert("[UNIMOD:172]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 1), ("O", 2), ("S", 1)])); // NBS
    composition.insert("[UNIMOD:176]".to_string(), HashMap::from([("H", 22), ("C", 15), ("O", 1)])); // BHT
    composition.insert("[UNIMOD:178]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 1), ("O", -1), ("S", 1)])); // DAET
    composition.insert("[UNIMOD:184]".to_string(), HashMap::from([("C", -9), ("13C", 9)])); // Label:13C(9)
    composition.insert("[UNIMOD:185]".to_string(), HashMap::from([("H", 1), ("C", -9), ("13C", 9), ("O", 3), ("P", 1)])); // Label:13C(9)+Phospho
    composition.insert("[UNIMOD:186]".to_string(), HashMap::from([("H", 4), ("C", 8), ("O", 2)])); // HPG
    composition.insert("[UNIMOD:187]".to_string(), HashMap::from([("H", 10), ("C", 16), ("O", 5)])); // 2HPG
    composition.insert("[UNIMOD:188]".to_string(), HashMap::from([("C", -6), ("13C", 6)])); // Label:13C(6)
    composition.insert("[UNIMOD:193]".to_string(), HashMap::from([("O", -2), ("18O", 2)])); // Label:18O(2)
    composition.insert("[UNIMOD:194]".to_string(), HashMap::from([("H", 6), ("C", 10), ("N", 2), ("O", 1)])); // AccQTag
    composition.insert("[UNIMOD:195]".to_string(), HashMap::from([("H", 19), ("C", 9), ("N", 2), ("O", 1)])); // QAT
    composition.insert("[UNIMOD:196]".to_string(), HashMap::from([("H", 16), ("2H", 3), ("C", 9), ("N", 2), ("O", 1)])); // QAT:2H(3)
    composition.insert("[UNIMOD:197]".to_string(), HashMap::from([("H", 20), ("C", 10), ("N", 2), ("O", 1)])); // EQAT
    composition.insert("[UNIMOD:198]".to_string(), HashMap::from([("H", 15), ("2H", 5), ("C", 10), ("N", 2), ("O", 1)])); // EQAT:2H(5)
    composition.insert("[UNIMOD:199]".to_string(), HashMap::from([("2H", 4), ("C", 2)])); // Dimethyl:2H(4)
    composition.insert("[UNIMOD:200]".to_string(), HashMap::from([("H", 4), ("C", 2), ("O", -1), ("S", 2)])); // Ethanedithiol
    composition.insert("[UNIMOD:205]".to_string(), HashMap::from([("H", 6), ("C", 6), ("O", 1)])); // Delta:H(6)C(6)O(1)
    composition.insert("[UNIMOD:206]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 1)])); // Delta:H(4)C(3)O(1)
    composition.insert("[UNIMOD:207]".to_string(), HashMap::from([("H", 2), ("C", 3)])); // Delta:H(2)C(3)
    composition.insert("[UNIMOD:208]".to_string(), HashMap::from([("H", 4), ("C", 6)])); // Delta:H(4)C(6)
    composition.insert("[UNIMOD:209]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 2)])); // Delta:H(8)C(6)O(2)
    composition.insert("[UNIMOD:211]".to_string(), HashMap::from([("H", 7), ("C", 4), ("N", 1), ("O", 1)])); // NEIAA
    composition.insert("[UNIMOD:212]".to_string(), HashMap::from([("H", 2), ("2H", 5), ("C", 4), ("N", 1), ("O", 1)])); // NEIAA:2H(5)
    composition.insert("[UNIMOD:213]".to_string(), HashMap::from([("H", 21), ("C", 15), ("N", 5), ("O", 13), ("P", 2)])); // ADP-Ribosyl
    composition.insert("[UNIMOD:214]".to_string(), HashMap::from([("H", 12), ("C", 4), ("13C", 3), ("N", 1), ("15N", 1), ("O", 1)])); // iTRAQ4plex
    composition.insert("[UNIMOD:243]".to_string(), HashMap::from([("H", 13), ("C", 12), ("N", 2), ("O", 2), ("Br", 1)])); // IGBP
    composition.insert("[UNIMOD:253]".to_string(), HashMap::from([("H", 6), ("C", 4), ("O", 1)])); // Crotonaldehyde
    composition.insert("[UNIMOD:254]".to_string(), HashMap::from([("H", 2), ("C", 2)])); // Delta:H(2)C(2)
    composition.insert("[UNIMOD:255]".to_string(), HashMap::from([("H", 4), ("C", 2)])); // Delta:H(4)C(2)
    composition.insert("[UNIMOD:256]".to_string(), HashMap::from([("H", 4), ("C", 3)])); // Delta:H(4)C(3)
    composition.insert("[UNIMOD:258]".to_string(), HashMap::from([("O", -1), ("18O", 1)])); // Label:18O(1)
    composition.insert("[UNIMOD:259]".to_string(), HashMap::from([("C", -6), ("13C", 6), ("N", -2), ("15N", 2)])); // Label:13C(6)15N(2)
    composition.insert("[UNIMOD:260]".to_string(), HashMap::from([("H", 1), ("O", 2), ("P", 1), ("S", 1)])); // Thiophospho
    composition.insert("[UNIMOD:261]".to_string(), HashMap::from([("H", 5), ("C", 7), ("N", 1), ("O", 3), ("S", 2)])); // SPITC
    composition.insert("[UNIMOD:262]".to_string(), HashMap::from([("H", -3), ("2H", 3)])); // Label:2H(3)
    composition.insert("[UNIMOD:264]".to_string(), HashMap::from([("H", 7), ("C", 7), ("N", 1), ("O", -1), ("S", 1)])); // PET
    composition.insert("[UNIMOD:267]".to_string(), HashMap::from([("C", -6), ("13C", 6), ("N", -4), ("15N", 4)])); // Label:13C(6)15N(4)
    composition.insert("[UNIMOD:268]".to_string(), HashMap::from([("C", -5), ("13C", 5), ("N", -1), ("15N", 1)])); // Label:13C(5)15N(1)
    composition.insert("[UNIMOD:269]".to_string(), HashMap::from([("C", -9), ("13C", 9), ("N", -1), ("15N", 1)])); // Label:13C(9)15N(1)
    composition.insert("[UNIMOD:270]".to_string(), HashMap::from([("H", 22), ("C", 19), ("O", 7)])); // Cytopiloyne
    composition.insert("[UNIMOD:271]".to_string(), HashMap::from([("H", 24), ("C", 19), ("O", 8)])); // Cytopiloyne+water
    composition.insert("[UNIMOD:272]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 4), ("S", 1)])); // CAF
    composition.insert("[UNIMOD:275]".to_string(), HashMap::from([("H", -1), ("N", 1), ("O", 1)])); // Nitrosyl
    composition.insert("[UNIMOD:276]".to_string(), HashMap::from([("H", 9), ("C", 8), ("N", 1), ("O", 2), ("S", 1)])); // AEBS
    composition.insert("[UNIMOD:278]".to_string(), HashMap::from([("H", 4), ("C", 2), ("O", 1)])); // Ethanolyl
    composition.insert("[UNIMOD:280]".to_string(), HashMap::from([("H", 4), ("C", 2)])); // Ethyl
    composition.insert("[UNIMOD:281]".to_string(), HashMap::from([("H", 34), ("C", 21), ("N", 7), ("O", 16), ("P", 3), ("S", 1)])); // CoenzymeA
    composition.insert("[UNIMOD:284]".to_string(), HashMap::from([("2H", 2), ("C", 1)])); // Methyl:2H(2)
    composition.insert("[UNIMOD:285]".to_string(), HashMap::from([("H", 5), ("C", 6), ("N", 1), ("O", 2), ("S", 1)])); // SulfanilicAcid
    composition.insert("[UNIMOD:286]".to_string(), HashMap::from([("H", 5), ("13C", 6), ("N", 1), ("O", 2), ("S", 1)])); // SulfanilicAcid:13C(6)
    composition.insert("[UNIMOD:288]".to_string(), HashMap::from([("H", -2), ("O", 1)])); // Trp->Oxolactone
    composition.insert("[UNIMOD:289]".to_string(), HashMap::from([("H", 28), ("C", 16), ("N", 4), ("O", 3), ("S", 1)])); // Biotin-PEO-Amine
    composition.insert("[UNIMOD:290]".to_string(), HashMap::from([("H", 32), ("C", 19), ("N", 4), ("O", 3), ("S", 2)])); // Biotin-HPDP
    composition.insert("[UNIMOD:291]".to_string(), HashMap::from([("Hg", 1)])); // Delta:Hg(1)
    composition.insert("[UNIMOD:292]".to_string(), HashMap::from([("H", 11), ("C", 9), ("N", 2), ("O", 9), ("P", 1)])); // IodoU-AMP
    composition.insert("[UNIMOD:293]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 2), ("S", 1)])); // CAMthiopropanoyl
    composition.insert("[UNIMOD:294]".to_string(), HashMap::from([("H", 22), ("C", 14), ("N", 4), ("O", 3), ("S", 1)])); // IED-Biotin
    composition.insert("[UNIMOD:295]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 4)])); // dHex
    composition.insert("[UNIMOD:298]".to_string(), HashMap::from([("H", -1), ("2H", 3), ("C", 1)])); // Methyl:2H(3)
    composition.insert("[UNIMOD:299]".to_string(), HashMap::from([("C", 1), ("O", 2)])); // Carboxy
    composition.insert("[UNIMOD:301]".to_string(), HashMap::from([("H", 10), ("C", 10), ("N", 2), ("O", 2)])); // Bromobimane
    composition.insert("[UNIMOD:302]".to_string(), HashMap::from([("H", 6), ("C", 11), ("O", 2)])); // Menadione
    composition.insert("[UNIMOD:303]".to_string(), HashMap::from([("H", 4), ("C", 2), ("O", 1), ("S", 1)])); // DeStreak
    composition.insert("[UNIMOD:305]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 39)])); // dHex(1)Hex(3)HexNAc(4)
    composition.insert("[UNIMOD:307]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 44)])); // dHex(1)Hex(4)HexNAc(4)
    composition.insert("[UNIMOD:308]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 49)])); // dHex(1)Hex(5)HexNAc(4)
    composition.insert("[UNIMOD:309]".to_string(), HashMap::from([("H", 82), ("C", 50), ("N", 4), ("O", 35)])); // Hex(3)HexNAc(4)
    composition.insert("[UNIMOD:310]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 40)])); // Hex(4)HexNAc(4)
    composition.insert("[UNIMOD:311]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 45)])); // Hex(5)HexNAc(4)
    composition.insert("[UNIMOD:312]".to_string(), HashMap::from([("H", 5), ("C", 3), ("N", 1), ("O", 2), ("S", 1)])); // Cysteinyl
    composition.insert("[UNIMOD:313]".to_string(), HashMap::from([("H", -12), ("C", -6), ("N", -2), ("O", -1)])); // Lys-loss
    composition.insert("[UNIMOD:314]".to_string(), HashMap::from([("H", 5), ("C", 5), ("N", 1), ("O", 2)])); // Nmethylmaleimide
    composition.insert("[UNIMOD:316]".to_string(), HashMap::from([("H", 6), ("C", 6)])); // DimethylpyrroleAdduct
    composition.insert("[UNIMOD:318]".to_string(), HashMap::from([("H", 2), ("C", 5)])); // Delta:H(2)C(5)
    composition.insert("[UNIMOD:319]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 1)])); // Delta:H(2)C(3)O(1)
    composition.insert("[UNIMOD:320]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 3)])); // Nethylmaleimide+water
    composition.insert("[UNIMOD:323]".to_string(), HashMap::from([("H", 30), ("C", 31), ("N", 4), ("O", 6), ("S", 1), ("I", 1)])); // Xlink:B10621
    composition.insert("[UNIMOD:324]".to_string(), HashMap::from([("H", 5), ("C", 3), ("N", 1), ("S", 1)])); // Xlink:DTBP[87]
    composition.insert("[UNIMOD:325]".to_string(), HashMap::from([("H", 49), ("C", 27), ("N", 4), ("O", 5), ("P", 1), ("S", 1)])); // FP-Biotin
    composition.insert("[UNIMOD:327]".to_string(), HashMap::from([("H", 4), ("C", 2), ("O", -1), ("S", 1)])); // Delta:H(4)C(2)O(-1)S(1)
    composition.insert("[UNIMOD:329]".to_string(), HashMap::from([("H", -1), ("2H", 3), ("13C", 1)])); // Methyl:2H(3)13C(1)
    composition.insert("[UNIMOD:330]".to_string(), HashMap::from([("H", -2), ("2H", 6), ("13C", 2)])); // Dimethyl:2H(6)13C(2)
    composition.insert("[UNIMOD:332]".to_string(), HashMap::from([("H", 34), ("C", 19), ("N", 4), ("O", 5), ("P", 1), ("S", 3)])); // Thiophos-S-S-biotin
    composition.insert("[UNIMOD:333]".to_string(), HashMap::from([("H", 34), ("C", 19), ("N", 3), ("O", 5), ("P", 1), ("S", 1)])); // Can-FP-biotin
    composition.insert("[UNIMOD:335]".to_string(), HashMap::from([("H", 18), ("C", 9), ("O", 2)])); // HNE+Delta:H(2)
    composition.insert("[UNIMOD:337]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 1), ("O", -1)])); // Methylamine
    composition.insert("[UNIMOD:340]".to_string(), HashMap::from([("H", -1), ("Br", 1)])); // Bromo
    composition.insert("[UNIMOD:342]".to_string(), HashMap::from([("H", 1), ("N", 1)])); // Amino
    composition.insert("[UNIMOD:343]".to_string(), HashMap::from([("H", 13), ("C", 9), ("N", 1), ("O", 2), ("S", 1)])); // Argbiotinhydrazide
    composition.insert("[UNIMOD:344]".to_string(), HashMap::from([("H", -5), ("C", -1), ("N", -3), ("O", 1)])); // Arg->GluSA
    composition.insert("[UNIMOD:345]".to_string(), HashMap::from([("O", 3)])); // Trioxidation
    composition.insert("[UNIMOD:348]".to_string(), HashMap::from([("H", -1), ("C", -2), ("N", -1), ("O", 1)])); // His->Asn
    composition.insert("[UNIMOD:349]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", -2), ("O", 2)])); // His->Asp
    composition.insert("[UNIMOD:350]".to_string(), HashMap::from([("C", -1), ("O", 2)])); // Trp->Hydroxykynurenin
    composition.insert("[UNIMOD:351]".to_string(), HashMap::from([("C", -1), ("O", 1)])); // Trp->Kynurenin
    composition.insert("[UNIMOD:352]".to_string(), HashMap::from([("H", -3), ("N", -1), ("O", 1)])); // Lys->Allysine
    composition.insert("[UNIMOD:353]".to_string(), HashMap::from([("H", 15), ("C", 10), ("N", 3), ("O", 2), ("S", 1)])); // Lysbiotinhydrazide
    composition.insert("[UNIMOD:354]".to_string(), HashMap::from([("H", -1), ("N", 1), ("O", 2)])); // Nitro
    composition.insert("[UNIMOD:357]".to_string(), HashMap::from([("H", 18), ("C", 10), ("N", 4), ("O", 2), ("S", 1)])); // probiotinhydrazide
    composition.insert("[UNIMOD:359]".to_string(), HashMap::from([("H", -2), ("O", 1)])); // Pro->pyro-Glu
    composition.insert("[UNIMOD:360]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", -1)])); // Pro->Pyrrolidinone
    composition.insert("[UNIMOD:361]".to_string(), HashMap::from([("H", 16), ("C", 10), ("N", 4), ("O", 1), ("S", 1)])); // Thrbiotinhydrazide
    composition.insert("[UNIMOD:362]".to_string(), HashMap::from([("H", 13), ("C", 6), ("O", 3), ("P", 1)])); // Diisopropylphosphate
    composition.insert("[UNIMOD:363]".to_string(), HashMap::from([("H", 7), ("C", 3), ("O", 3), ("P", 1)])); // Isopropylphospho
    composition.insert("[UNIMOD:364]".to_string(), HashMap::from([("H", 3), ("13C", 6), ("N", 1), ("O", 1)])); // ICPL:13C(6)
    composition.insert("[UNIMOD:365]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 1), ("O", 1)])); // ICPL
    composition.insert("[UNIMOD:366]".to_string(), HashMap::from([("H", -1), ("N", -1), ("18O", 1)])); // Deamidated:18O(1)
    composition.insert("[UNIMOD:368]".to_string(), HashMap::from([("H", -2), ("S", -1)])); // Cys->Dha
    composition.insert("[UNIMOD:369]".to_string(), HashMap::from([("C", -1), ("O", -1)])); // Pro->Pyrrolidone
    composition.insert("[UNIMOD:371]".to_string(), HashMap::from([("H", 6), ("C", 4), ("O", 2)])); // HMVK
    composition.insert("[UNIMOD:372]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", -2)])); // Arg->Orn
    composition.insert("[UNIMOD:374]".to_string(), HashMap::from([("H", -1)])); // Dehydro
    composition.insert("[UNIMOD:375]".to_string(), HashMap::from([("H", 14), ("C", 7), ("N", 2), ("O", 1)])); // Diphthamide
    composition.insert("[UNIMOD:376]".to_string(), HashMap::from([("H", 24), ("C", 15), ("O", 1)])); // Hydroxyfarnesyl
    composition.insert("[UNIMOD:377]".to_string(), HashMap::from([("H", 68), ("C", 37), ("O", 4)])); // Diacylglycerol
    composition.insert("[UNIMOD:378]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 2)])); // Carboxyethyl
    composition.insert("[UNIMOD:379]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 1), ("O", 1)])); // Hypusine
    composition.insert("[UNIMOD:380]".to_string(), HashMap::from([("H", 26), ("C", 20)])); // Retinylidene
    composition.insert("[UNIMOD:381]".to_string(), HashMap::from([("H", -3), ("N", -1), ("O", 2)])); // Lys->AminoadipicAcid
    composition.insert("[UNIMOD:382]".to_string(), HashMap::from([("H", -3), ("N", -1), ("O", 1), ("S", -1)])); // Cys->PyruvicAcid
    composition.insert("[UNIMOD:385]".to_string(), HashMap::from([("H", -3), ("N", -1)])); // Ammonia-loss
    composition.insert("[UNIMOD:387]".to_string(), HashMap::from([("H", 38), ("C", 33), ("N", 4), ("O", 6)])); // Phycocyanobilin
    composition.insert("[UNIMOD:388]".to_string(), HashMap::from([("H", 40), ("C", 33), ("N", 4), ("O", 6)])); // Phycoerythrobilin
    composition.insert("[UNIMOD:389]".to_string(), HashMap::from([("H", 36), ("C", 33), ("N", 4), ("O", 6)])); // Phytochromobilin
    composition.insert("[UNIMOD:390]".to_string(), HashMap::from([("H", 32), ("C", 34), ("N", 4), ("O", 4), ("Fe", 1)])); // Heme
    composition.insert("[UNIMOD:391]".to_string(), HashMap::from([("H", 11), ("C", 10), ("N", 5), ("O", 8), ("P", 1), ("S", 2), ("Mo", 1)])); // Molybdopterin
    composition.insert("[UNIMOD:392]".to_string(), HashMap::from([("H", -2), ("O", 2)])); // Quinone
    composition.insert("[UNIMOD:393]".to_string(), HashMap::from([("H", 20), ("C", 12), ("O", 11)])); // Glucosylgalactosyl
    composition.insert("[UNIMOD:394]".to_string(), HashMap::from([("H", 6), ("C", 2), ("N", 1), ("O", 3), ("P", 1)])); // GPIanchor
    composition.insert("[UNIMOD:395]".to_string(), HashMap::from([("H", 42), ("C", 26), ("N", 7), ("O", 19), ("P", 3), ("S", 1)])); // PhosphoribosyldephosphoCoA
    composition.insert("[UNIMOD:396]".to_string(), HashMap::from([("H", 12), ("C", 5), ("N", 1), ("O", 5), ("P", 1)])); // GlycerylPE
    composition.insert("[UNIMOD:397]".to_string(), HashMap::from([("H", 1), ("C", 6), ("O", 1), ("I", 3)])); // Triiodothyronine
    composition.insert("[UNIMOD:398]".to_string(), HashMap::from([("C", 6), ("O", 1), ("I", 4)])); // Thyroxine
    composition.insert("[UNIMOD:400]".to_string(), HashMap::from([("H", -6), ("C", -6), ("O", -1)])); // Tyr->Dha
    composition.insert("[UNIMOD:401]".to_string(), HashMap::from([("H", -2)])); // Didehydro
    composition.insert("[UNIMOD:402]".to_string(), HashMap::from([("H", -2), ("O", 1), ("S", -1)])); // Cys->Oxoalanine
    composition.insert("[UNIMOD:403]".to_string(), HashMap::from([("H", -1), ("N", -1)])); // Ser->LacticAcid
    composition.insert("[UNIMOD:405]".to_string(), HashMap::from([("H", 12), ("C", 10), ("N", 5), ("O", 6), ("P", 1)])); // Phosphoadenosine
    composition.insert("[UNIMOD:407]".to_string(), HashMap::from([("H", 6), ("C", 9), ("O", 2)])); // Hydroxycinnamyl
    composition.insert("[UNIMOD:408]".to_string(), HashMap::from([("H", 8), ("C", 5), ("O", 5)])); // Glycosyl
    composition.insert("[UNIMOD:409]".to_string(), HashMap::from([("H", 19), ("C", 17), ("N", 4), ("O", 9), ("P", 1)])); // FMNH
    composition.insert("[UNIMOD:410]".to_string(), HashMap::from([("H", 86), ("C", 43), ("O", 2)])); // Archaeol
    composition.insert("[UNIMOD:411]".to_string(), HashMap::from([("H", 5), ("C", 7), ("N", 1), ("O", 1)])); // Phenylisocyanate
    composition.insert("[UNIMOD:412]".to_string(), HashMap::from([("2H", 5), ("C", 7), ("N", 1), ("O", 1)])); // Phenylisocyanate:2H(5)
    composition.insert("[UNIMOD:413]".to_string(), HashMap::from([("H", 12), ("C", 10), ("N", 5), ("O", 7), ("P", 1)])); // Phosphoguanosine
    composition.insert("[UNIMOD:414]".to_string(), HashMap::from([("H", 2), ("C", 1), ("O", 1)])); // Hydroxymethyl
    composition.insert("[UNIMOD:415]".to_string(), HashMap::from([("H", 47), ("C", 40), ("N", 20), ("O", 26), ("P", 4), ("S", 3), ("Se", 1), ("Mo", 1)])); // MolybdopterinGD+Delta:S(-1)Se(1)
    composition.insert("[UNIMOD:416]".to_string(), HashMap::from([("H", 22), ("C", 20), ("N", 2), ("O", 8)])); // Dipyrrolylmethanemethyl
    composition.insert("[UNIMOD:417]".to_string(), HashMap::from([("H", 11), ("C", 9), ("N", 2), ("O", 8), ("P", 1)])); // PhosphoUridine
    composition.insert("[UNIMOD:419]".to_string(), HashMap::from([("H", 7), ("C", 3), ("O", 5), ("P", 1)])); // Glycerophospho
    composition.insert("[UNIMOD:420]".to_string(), HashMap::from([("O", -1), ("S", 1)])); // Carboxy->Thiocarboxy
    composition.insert("[UNIMOD:421]".to_string(), HashMap::from([("S", 1)])); // Sulfide
    composition.insert("[UNIMOD:422]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 2)])); // PyruvicAcidIminyl
    composition.insert("[UNIMOD:423]".to_string(), HashMap::from([("Se", 1)])); // Delta:Se(1)
    composition.insert("[UNIMOD:424]".to_string(), HashMap::from([("H", 47), ("C", 40), ("N", 20), ("O", 26), ("P", 4), ("S", 4), ("Mo", 1)])); // MolybdopterinGD
    composition.insert("[UNIMOD:425]".to_string(), HashMap::from([("O", 2)])); // Dioxidation
    composition.insert("[UNIMOD:426]".to_string(), HashMap::from([("H", 14), ("C", 8), ("O", 1)])); // Octanoyl
    composition.insert("[UNIMOD:428]".to_string(), HashMap::from([("H", 14), ("C", 8), ("N", 1), ("O", 8), ("P", 1)])); // PhosphoHexNAc
    composition.insert("[UNIMOD:429]".to_string(), HashMap::from([("H", 11), ("C", 6), ("O", 8), ("P", 1)])); // PhosphoHex
    composition.insert("[UNIMOD:431]".to_string(), HashMap::from([("H", 28), ("C", 16), ("O", 1)])); // Palmitoleyl
    composition.insert("[UNIMOD:432]".to_string(), HashMap::from([("H", 44), ("C", 27)])); // Cholesterol
    composition.insert("[UNIMOD:433]".to_string(), HashMap::from([("H", 24), ("C", 20)])); // Didehydroretinylidene
    composition.insert("[UNIMOD:434]".to_string(), HashMap::from([("H", 26), ("C", 17), ("O", 4)])); // CHDH
    composition.insert("[UNIMOD:435]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 1), ("O", 1)])); // Methylpyrroline
    composition.insert("[UNIMOD:436]".to_string(), HashMap::from([("H", 30), ("C", 34), ("N", 4), ("O", 4), ("Fe", 1)])); // Hydroxyheme
    composition.insert("[UNIMOD:437]".to_string(), HashMap::from([("H", 19), ("C", 13), ("N", 6), ("O", 6), ("P", 1)])); // MicrocinC7
    composition.insert("[UNIMOD:438]".to_string(), HashMap::from([("H", -1), ("C", 1), ("N", 1)])); // Cyano
    composition.insert("[UNIMOD:439]".to_string(), HashMap::from([("H", -1), ("C", 5), ("N", 2), ("O", 5), ("S", 2), ("Fe", 2)])); // Diironsubcluster
    composition.insert("[UNIMOD:440]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 2)])); // Amidino
    composition.insert("[UNIMOD:442]".to_string(), HashMap::from([("H", 19), ("C", 17), ("N", 4), ("O", 8), ("P", 1)])); // FMN
    composition.insert("[UNIMOD:443]".to_string(), HashMap::from([("H", 21), ("C", 17), ("N", 4), ("O", 9), ("P", 1)])); // FMNC
    composition.insert("[UNIMOD:444]".to_string(), HashMap::from([("H", 24), ("C", 19), ("N", 8), ("O", 15), ("P", 2), ("S", 3), ("Mo", 1), ("Cu", 1)])); // CuSMo
    composition.insert("[UNIMOD:445]".to_string(), HashMap::from([("H", 7), ("C", 3), ("O", 1)])); // Hydroxytrimethyl
    composition.insert("[UNIMOD:447]".to_string(), HashMap::from([("O", -1)])); // Deoxy
    composition.insert("[UNIMOD:448]".to_string(), HashMap::from([("H", 37), ("C", 36), ("N", 3), ("O", 20)])); // Microcin
    composition.insert("[UNIMOD:449]".to_string(), HashMap::from([("H", 18), ("C", 10), ("O", 1)])); // Decanoyl
    composition.insert("[UNIMOD:450]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 3)])); // Glu
    composition.insert("[UNIMOD:451]".to_string(), HashMap::from([("H", 14), ("C", 10), ("N", 2), ("O", 6)])); // GluGlu
    composition.insert("[UNIMOD:452]".to_string(), HashMap::from([("H", 21), ("C", 15), ("N", 3), ("O", 9)])); // GluGluGlu
    composition.insert("[UNIMOD:453]".to_string(), HashMap::from([("H", 28), ("C", 20), ("N", 4), ("O", 12)])); // GluGluGluGlu
    composition.insert("[UNIMOD:454]".to_string(), HashMap::from([("H", 11), ("C", 6), ("N", 1), ("O", 4)])); // HexN
    composition.insert("[UNIMOD:455]".to_string(), HashMap::from([("H", 14), ("C", 8), ("N", 2), ("O", 1)])); // Xlink:DMP[154]
    composition.insert("[UNIMOD:457]".to_string(), HashMap::from([("H", 5), ("C", 13), ("N", 1)])); // NDA
    composition.insert("[UNIMOD:464]".to_string(), HashMap::from([("H", 5), ("C", 1), ("13C", 6), ("N", 1), ("O", 3), ("S", 2)])); // SPITC:13C(6)
    composition.insert("[UNIMOD:472]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 1), ("O", -1), ("S", 1)])); // AEC-MAEC
    composition.insert("[UNIMOD:476]".to_string(), HashMap::from([("H", 14), ("C", 7), ("N", 1), ("O", 1)])); // TMAB
    composition.insert("[UNIMOD:477]".to_string(), HashMap::from([("H", 5), ("2H", 9), ("C", 7), ("N", 1), ("O", 1)])); // TMAB:2H(9)
    composition.insert("[UNIMOD:478]".to_string(), HashMap::from([("H", 15), ("C", 21), ("N", 3), ("O", 5), ("S", 1)])); // FTC
    composition.insert("[UNIMOD:481]".to_string(), HashMap::from([("H", -4), ("2H", 4)])); // Label:2H(4)
    composition.insert("[UNIMOD:488]".to_string(), HashMap::from([("H", 8), ("C", 8), ("N", 1)])); // DHP
    composition.insert("[UNIMOD:490]".to_string(), HashMap::from([("H", 12), ("C", 7), ("O", 6)])); // Hep
    composition.insert("[UNIMOD:493]".to_string(), HashMap::from([("H", 24), ("C", 21), ("O", 4)])); // BADGE
    composition.insert("[UNIMOD:494]".to_string(), HashMap::from([("H", 44), ("C", 37), ("N", 4), ("O", 6), ("S", 1)])); // CyDye-Cy3
    composition.insert("[UNIMOD:495]".to_string(), HashMap::from([("H", 44), ("C", 38), ("N", 4), ("O", 6), ("S", 1)])); // CyDye-Cy5
    composition.insert("[UNIMOD:498]".to_string(), HashMap::from([("H", 22), ("C", 15), ("O", 2)])); // BHTOH
    composition.insert("[UNIMOD:499]".to_string(), HashMap::from([("H", 13), ("C", 10), ("13C", 2), ("N", 2), ("O", 2), ("Br", 1)])); // IGBP:13C(2)
    composition.insert("[UNIMOD:500]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 3)])); // Nmethylmaleimide+water
    composition.insert("[UNIMOD:501]".to_string(), HashMap::from([("H", 6), ("C", 7), ("N", 2), ("O", 1)])); // PyMIC
    composition.insert("[UNIMOD:503]".to_string(), HashMap::from([("H", 28), ("C", 20), ("O", 4)])); // LG-lactam-K
    composition.insert("[UNIMOD:504]".to_string(), HashMap::from([("H", 28), ("C", 20), ("O", 5)])); // LG-Hlactam-K
    composition.insert("[UNIMOD:505]".to_string(), HashMap::from([("H", 26), ("C", 19), ("N", -2), ("O", 4)])); // LG-lactam-R
    composition.insert("[UNIMOD:506]".to_string(), HashMap::from([("H", 26), ("C", 19), ("N", -2), ("O", 5)])); // LG-Hlactam-R
    composition.insert("[UNIMOD:510]".to_string(), HashMap::from([("2H", 4), ("13C", 2)])); // Dimethyl:2H(4)13C(2)
    composition.insert("[UNIMOD:512]".to_string(), HashMap::from([("H", 20), ("C", 12), ("O", 10)])); // Hex(2)
    composition.insert("[UNIMOD:513]".to_string(), HashMap::from([("H", 29), ("C", 14), ("N", 1), ("O", 1)])); // C8-QAT
    composition.insert("[UNIMOD:514]".to_string(), HashMap::from([("H", 14), ("C", 9), ("N", 1), ("O", 4), ("S", 1)])); // PropylNAGthiazoline
    composition.insert("[UNIMOD:515]".to_string(), HashMap::from([("H", 13), ("C", 24), ("N", 1), ("O", 7)])); // FNEM
    composition.insert("[UNIMOD:518]".to_string(), HashMap::from([("H", 8), ("C", 4)])); // Diethyl
    composition.insert("[UNIMOD:519]".to_string(), HashMap::from([("H", 22), ("C", 32), ("N", 2), ("O", 6), ("S", 2)])); // BisANS
    composition.insert("[UNIMOD:520]".to_string(), HashMap::from([("H", 8), ("C", 5)])); // Piperidine
    composition.insert("[UNIMOD:522]".to_string(), HashMap::from([("H", 35), ("C", 23), ("N", 5), ("O", 7), ("S", 1)])); // Maleimide-PEO2-Biotin
    composition.insert("[UNIMOD:523]".to_string(), HashMap::from([("H", 36), ("C", 22), ("N", 4), ("O", 4), ("S", 1)])); // Sulfo-NHS-LC-LC-Biotin
    composition.insert("[UNIMOD:525]".to_string(), HashMap::from([("H", 12), ("C", 6), ("13C", 1), ("N", 2), ("O", 1)])); // CLIP_TRAQ_2
    composition.insert("[UNIMOD:526]".to_string(), HashMap::from([("H", -4), ("C", -1), ("S", -1)])); // Dethiomethyl
    composition.insert("[UNIMOD:528]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", -1), ("O", 1)])); // Methyl+Deamidated
    composition.insert("[UNIMOD:529]".to_string(), HashMap::from([("H", 5), ("C", 2)])); // Delta:H(5)C(2)
    composition.insert("[UNIMOD:530]".to_string(), HashMap::from([("H", -1), ("K", 1)])); // Cation:K
    composition.insert("[UNIMOD:531]".to_string(), HashMap::from([("H", -1), ("Cu", 1)])); // Cation:Cu[I]
    composition.insert("[UNIMOD:532]".to_string(), HashMap::from([("H", 12), ("C", 5), ("13C", 2), ("N", 2), ("18O", 1)])); // iTRAQ4plex114
    composition.insert("[UNIMOD:533]".to_string(), HashMap::from([("H", 12), ("C", 6), ("13C", 1), ("N", 1), ("15N", 1), ("18O", 1)])); // iTRAQ4plex115
    composition.insert("[UNIMOD:534]".to_string(), HashMap::from([("H", -2), ("Br", 2)])); // Dibromo
    composition.insert("[UNIMOD:535]".to_string(), HashMap::from([("H", 29), ("C", 16), ("N", 7), ("O", 4)])); // LRGG
    composition.insert("[UNIMOD:536]".to_string(), HashMap::from([("H", 20), ("C", 11), ("13C", 1), ("N", 3), ("O", 4)])); // CLIP_TRAQ_3
    composition.insert("[UNIMOD:537]".to_string(), HashMap::from([("H", 15), ("C", 9), ("13C", 1), ("N", 2), ("O", 5)])); // CLIP_TRAQ_4
    composition.insert("[UNIMOD:538]".to_string(), HashMap::from([("H", 54), ("C", 35), ("N", 4), ("O", 4), ("S", 1)])); // Biotin:Cayman-10141
    composition.insert("[UNIMOD:539]".to_string(), HashMap::from([("H", 60), ("C", 36), ("N", 4), ("O", 5), ("S", 1)])); // Biotin:Cayman-10013
    composition.insert("[UNIMOD:540]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 1), ("S", 0)])); // Ala->Ser
    composition.insert("[UNIMOD:541]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 1), ("S", 0)])); // Ala->Thr
    composition.insert("[UNIMOD:542]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 0), ("O", 2), ("S", 0)])); // Ala->Asp
    composition.insert("[UNIMOD:543]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 0), ("S", 0)])); // Ala->Pro
    composition.insert("[UNIMOD:544]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 0), ("S", 0)])); // Ala->Gly
    composition.insert("[UNIMOD:545]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 2), ("S", 0)])); // Ala->Glu
    composition.insert("[UNIMOD:546]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", 0), ("S", 0)])); // Ala->Val
    composition.insert("[UNIMOD:547]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", 0), ("S", -1)])); // Cys->Phe
    composition.insert("[UNIMOD:548]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 1), ("S", -1)])); // Cys->Ser
    composition.insert("[UNIMOD:549]".to_string(), HashMap::from([("H", 5), ("C", 8), ("N", 1), ("O", 0), ("S", -1)])); // Cys->Trp
    composition.insert("[UNIMOD:550]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", 1), ("S", -1)])); // Cys->Tyr
    composition.insert("[UNIMOD:551]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 3), ("O", 0), ("S", -1)])); // Cys->Arg
    composition.insert("[UNIMOD:552]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 0), ("S", -1)])); // Cys->Gly
    composition.insert("[UNIMOD:553]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", 0), ("O", -2), ("S", 0)])); // Asp->Ala
    composition.insert("[UNIMOD:554]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 2), ("O", -2), ("S", 0)])); // Asp->His
    composition.insert("[UNIMOD:555]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", 1), ("O", -1), ("S", 0)])); // Asp->Asn
    composition.insert("[UNIMOD:556]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", -2), ("S", 0)])); // Asp->Gly
    composition.insert("[UNIMOD:557]".to_string(), HashMap::from([("H", 4), ("C", 5), ("N", 0), ("O", -1), ("S", 0)])); // Asp->Tyr
    composition.insert("[UNIMOD:558]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 0), ("S", 0)])); // Asp->Glu
    composition.insert("[UNIMOD:559]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 0), ("O", -2), ("S", 0)])); // Asp->Val
    composition.insert("[UNIMOD:560]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", -2), ("S", 0)])); // Glu->Ala
    composition.insert("[UNIMOD:561]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", 1), ("O", -1), ("S", 0)])); // Glu->Gln
    composition.insert("[UNIMOD:562]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 0), ("S", 0)])); // Glu->Asp
    composition.insert("[UNIMOD:563]".to_string(), HashMap::from([("H", 5), ("C", 1), ("N", 1), ("O", -2), ("S", 0)])); // Glu->Lys
    composition.insert("[UNIMOD:564]".to_string(), HashMap::from([("H", -4), ("C", -3), ("N", 0), ("O", -2), ("S", 0)])); // Glu->Gly
    composition.insert("[UNIMOD:565]".to_string(), HashMap::from([("H", 2), ("C", 0), ("N", 0), ("O", -2), ("S", 0)])); // Glu->Val
    composition.insert("[UNIMOD:566]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", 1), ("S", 0)])); // Phe->Ser
    composition.insert("[UNIMOD:567]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", 0), ("S", 1)])); // Phe->Cys
    composition.insert("[UNIMOD:568]".to_string(), HashMap::from([("H", 2), ("C", -3)])); // Phe->Xle
    composition.insert("[UNIMOD:569]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 1), ("S", 0)])); // Phe->Tyr
    composition.insert("[UNIMOD:570]".to_string(), HashMap::from([("H", 0), ("C", -4), ("N", 0), ("O", 0), ("S", 0)])); // Phe->Val
    composition.insert("[UNIMOD:571]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 0), ("S", 0)])); // Gly->Ala
    composition.insert("[UNIMOD:572]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 1), ("S", 0)])); // Gly->Ser
    composition.insert("[UNIMOD:573]".to_string(), HashMap::from([("H", 7), ("C", 9), ("N", 1), ("O", 0), ("S", 0)])); // Gly->Trp
    composition.insert("[UNIMOD:574]".to_string(), HashMap::from([("H", 4), ("C", 3), ("N", 0), ("O", 2), ("S", 0)])); // Gly->Glu
    composition.insert("[UNIMOD:575]".to_string(), HashMap::from([("H", 6), ("C", 3), ("N", 0), ("O", 0), ("S", 0)])); // Gly->Val
    composition.insert("[UNIMOD:576]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 2), ("S", 0)])); // Gly->Asp
    composition.insert("[UNIMOD:577]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 0), ("S", 1)])); // Gly->Cys
    composition.insert("[UNIMOD:578]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 3), ("O", 0), ("S", 0)])); // Gly->Arg
    composition.insert("[UNIMOD:580]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", -2), ("O", 0), ("S", 0)])); // His->Pro
    composition.insert("[UNIMOD:581]".to_string(), HashMap::from([("H", 2), ("C", 3), ("N", -2), ("O", 1), ("S", 0)])); // His->Tyr
    composition.insert("[UNIMOD:582]".to_string(), HashMap::from([("H", 1), ("C", -1), ("N", -1), ("O", 1), ("S", 0)])); // His->Gln
    composition.insert("[UNIMOD:584]".to_string(), HashMap::from([("H", 5), ("C", 0), ("N", 1), ("O", 0), ("S", 0)])); // His->Arg
    composition.insert("[UNIMOD:585]".to_string(), HashMap::from([("H", 4), ("N", -2)])); // His->Xle
    composition.insert("[UNIMOD:588]".to_string(), HashMap::from([("H", -4), ("C", -2), ("O", 1)])); // Xle->Thr
    composition.insert("[UNIMOD:589]".to_string(), HashMap::from([("H", -5), ("C", -2), ("N", 1), ("O", 1)])); // Xle->Asn
    composition.insert("[UNIMOD:590]".to_string(), HashMap::from([("H", 1), ("N", 1)])); // Xle->Lys
    composition.insert("[UNIMOD:594]".to_string(), HashMap::from([("H", -5), ("C", -2), ("N", -1), ("O", 1), ("S", 0)])); // Lys->Thr
    composition.insert("[UNIMOD:595]".to_string(), HashMap::from([("H", -6), ("C", -2), ("N", 0), ("O", 1), ("S", 0)])); // Lys->Asn
    composition.insert("[UNIMOD:596]".to_string(), HashMap::from([("H", -5), ("C", -1), ("N", -1), ("O", 2), ("S", 0)])); // Lys->Glu
    composition.insert("[UNIMOD:597]".to_string(), HashMap::from([("H", -4), ("C", -1), ("N", 0), ("O", 1), ("S", 0)])); // Lys->Gln
    composition.insert("[UNIMOD:598]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -1), ("O", 0), ("S", 1)])); // Lys->Met
    composition.insert("[UNIMOD:599]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 2), ("O", 0), ("S", 0)])); // Lys->Arg
    composition.insert("[UNIMOD:600]".to_string(), HashMap::from([("H", -1), ("N", -1)])); // Lys->Xle
    composition.insert("[UNIMOD:601]".to_string(), HashMap::from([("H", -6), ("C", -3), ("O", 1)])); // Xle->Ser
    composition.insert("[UNIMOD:602]".to_string(), HashMap::from([("H", -2), ("C", 3)])); // Xle->Phe
    composition.insert("[UNIMOD:603]".to_string(), HashMap::from([("H", -1), ("C", 5), ("N", 1)])); // Xle->Trp
    composition.insert("[UNIMOD:604]".to_string(), HashMap::from([("H", -4), ("C", -1)])); // Xle->Pro
    composition.insert("[UNIMOD:605]".to_string(), HashMap::from([("H", -2), ("C", -1)])); // Xle->Val
    composition.insert("[UNIMOD:606]".to_string(), HashMap::from([("H", -4), ("N", 2)])); // Xle->His
    composition.insert("[UNIMOD:607]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", 1), ("O", 1)])); // Xle->Gln
    composition.insert("[UNIMOD:608]".to_string(), HashMap::from([("H", -2), ("C", -1), ("S", 1)])); // Xle->Met
    composition.insert("[UNIMOD:609]".to_string(), HashMap::from([("H", 1), ("N", 3)])); // Xle->Arg
    composition.insert("[UNIMOD:610]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 1), ("S", -1)])); // Met->Thr
    composition.insert("[UNIMOD:611]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 3), ("O", 0), ("S", -1)])); // Met->Arg
    composition.insert("[UNIMOD:613]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 1), ("O", 0), ("S", -1)])); // Met->Lys
    composition.insert("[UNIMOD:614]".to_string(), HashMap::from([("H", 2), ("C", 1), ("S", -1)])); // Met->Xle
    composition.insert("[UNIMOD:615]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 0), ("S", -1)])); // Met->Val
    composition.insert("[UNIMOD:616]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", -1), ("O", 0), ("S", 0)])); // Asn->Ser
    composition.insert("[UNIMOD:617]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", -1), ("O", 0), ("S", 0)])); // Asn->Thr
    composition.insert("[UNIMOD:618]".to_string(), HashMap::from([("H", 6), ("C", 2), ("N", 0), ("O", -1), ("S", 0)])); // Asn->Lys
    composition.insert("[UNIMOD:619]".to_string(), HashMap::from([("H", 3), ("C", 5), ("N", -1), ("O", 0), ("S", 0)])); // Asn->Tyr
    composition.insert("[UNIMOD:620]".to_string(), HashMap::from([("H", 1), ("C", 2), ("N", 1), ("O", -1), ("S", 0)])); // Asn->His
    composition.insert("[UNIMOD:621]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", -1), ("O", 1), ("S", 0)])); // Asn->Asp
    composition.insert("[UNIMOD:622]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", -1), ("O", -1)])); // Asn->Xle
    composition.insert("[UNIMOD:623]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", 1), ("S", 0)])); // Pro->Ser
    composition.insert("[UNIMOD:624]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", 0), ("S", 0)])); // Pro->Ala
    composition.insert("[UNIMOD:625]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 2), ("O", 0), ("S", 0)])); // Pro->His
    composition.insert("[UNIMOD:626]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", 1), ("O", 1), ("S", 0)])); // Pro->Gln
    composition.insert("[UNIMOD:627]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", 0), ("O", 1), ("S", 0)])); // Pro->Thr
    composition.insert("[UNIMOD:628]".to_string(), HashMap::from([("H", 5), ("C", 1), ("N", 3), ("O", 0), ("S", 0)])); // Pro->Arg
    composition.insert("[UNIMOD:629]".to_string(), HashMap::from([("H", 4), ("C", 1)])); // Pro->Xle
    composition.insert("[UNIMOD:630]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", -1), ("O", -1), ("S", 0)])); // Gln->Pro
    composition.insert("[UNIMOD:631]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 0), ("O", -1), ("S", 0)])); // Gln->Lys
    composition.insert("[UNIMOD:632]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", -1), ("O", 1), ("S", 0)])); // Gln->Glu
    composition.insert("[UNIMOD:633]".to_string(), HashMap::from([("H", -1), ("C", 1), ("N", 1), ("O", -1), ("S", 0)])); // Gln->His
    composition.insert("[UNIMOD:634]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 2), ("O", -1), ("S", 0)])); // Gln->Arg
    composition.insert("[UNIMOD:635]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", -1), ("O", -1)])); // Gln->Xle
    composition.insert("[UNIMOD:636]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -3), ("O", 1), ("S", 0)])); // Arg->Ser
    composition.insert("[UNIMOD:637]".to_string(), HashMap::from([("H", -2), ("C", 5), ("N", -2), ("O", 0), ("S", 0)])); // Arg->Trp
    composition.insert("[UNIMOD:638]".to_string(), HashMap::from([("H", -5), ("C", -2), ("N", -3), ("O", 1), ("S", 0)])); // Arg->Thr
    composition.insert("[UNIMOD:639]".to_string(), HashMap::from([("H", -5), ("C", -1), ("N", -3), ("O", 0), ("S", 0)])); // Arg->Pro
    composition.insert("[UNIMOD:640]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", -2), ("O", 0), ("S", 0)])); // Arg->Lys
    composition.insert("[UNIMOD:641]".to_string(), HashMap::from([("H", -5), ("C", 0), ("N", -1), ("O", 0), ("S", 0)])); // Arg->His
    composition.insert("[UNIMOD:642]".to_string(), HashMap::from([("H", -4), ("C", -1), ("N", -2), ("O", 1), ("S", 0)])); // Arg->Gln
    composition.insert("[UNIMOD:643]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -3), ("O", 0), ("S", 1)])); // Arg->Met
    composition.insert("[UNIMOD:644]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -3), ("O", 0), ("S", 1)])); // Arg->Cys
    composition.insert("[UNIMOD:645]".to_string(), HashMap::from([("H", -1), ("N", -3)])); // Arg->Xle
    composition.insert("[UNIMOD:646]".to_string(), HashMap::from([("H", -9), ("C", -4), ("N", -3), ("O", 0), ("S", 0)])); // Arg->Gly
    composition.insert("[UNIMOD:647]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", -1), ("S", 0)])); // Ser->Phe
    composition.insert("[UNIMOD:648]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", -1), ("S", 0)])); // Ser->Ala
    composition.insert("[UNIMOD:649]".to_string(), HashMap::from([("H", 5), ("C", 8), ("N", 1), ("O", -1), ("S", 0)])); // Ser->Trp
    composition.insert("[UNIMOD:650]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 0), ("S", 0)])); // Ser->Thr
    composition.insert("[UNIMOD:651]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", 1), ("O", 0), ("S", 0)])); // Ser->Asn
    composition.insert("[UNIMOD:652]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", -1), ("S", 0)])); // Ser->Pro
    composition.insert("[UNIMOD:653]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", 0), ("S", 0)])); // Ser->Tyr
    composition.insert("[UNIMOD:654]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", -1), ("S", 1)])); // Ser->Cys
    composition.insert("[UNIMOD:655]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 3), ("O", -1), ("S", 0)])); // Ser->Arg
    composition.insert("[UNIMOD:656]".to_string(), HashMap::from([("H", 6), ("C", 3), ("O", -1)])); // Ser->Xle
    composition.insert("[UNIMOD:657]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", -1), ("S", 0)])); // Ser->Gly
    composition.insert("[UNIMOD:658]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 0), ("S", 0)])); // Thr->Ser
    composition.insert("[UNIMOD:659]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", -1), ("S", 0)])); // Thr->Ala
    composition.insert("[UNIMOD:660]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", 1), ("O", 0), ("S", 0)])); // Thr->Asn
    composition.insert("[UNIMOD:661]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 1), ("O", -1), ("S", 0)])); // Thr->Lys
    composition.insert("[UNIMOD:662]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 0), ("O", -1), ("S", 0)])); // Thr->Pro
    composition.insert("[UNIMOD:663]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", -1), ("S", 1)])); // Thr->Met
    composition.insert("[UNIMOD:664]".to_string(), HashMap::from([("H", 4), ("C", 2), ("O", -1)])); // Thr->Xle
    composition.insert("[UNIMOD:665]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 3), ("O", -1), ("S", 0)])); // Thr->Arg
    composition.insert("[UNIMOD:666]".to_string(), HashMap::from([("H", 0), ("C", 4), ("N", 0), ("O", 0), ("S", 0)])); // Val->Phe
    composition.insert("[UNIMOD:667]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 0), ("S", 0)])); // Val->Ala
    composition.insert("[UNIMOD:668]".to_string(), HashMap::from([("H", -2), ("C", 0), ("N", 0), ("O", 2), ("S", 0)])); // Val->Glu
    composition.insert("[UNIMOD:669]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 0), ("S", 1)])); // Val->Met
    composition.insert("[UNIMOD:670]".to_string(), HashMap::from([("H", -4), ("C", -1), ("N", 0), ("O", 2), ("S", 0)])); // Val->Asp
    composition.insert("[UNIMOD:671]".to_string(), HashMap::from([("H", 2), ("C", 1)])); // Val->Xle
    composition.insert("[UNIMOD:672]".to_string(), HashMap::from([("H", -6), ("C", -3), ("N", 0), ("O", 0), ("S", 0)])); // Val->Gly
    composition.insert("[UNIMOD:673]".to_string(), HashMap::from([("H", -5), ("C", -8), ("N", -1), ("O", 1), ("S", 0)])); // Trp->Ser
    composition.insert("[UNIMOD:674]".to_string(), HashMap::from([("H", -5), ("C", -8), ("N", -1), ("O", 0), ("S", 1)])); // Trp->Cys
    composition.insert("[UNIMOD:675]".to_string(), HashMap::from([("H", 2), ("C", -5), ("N", 2), ("O", 0), ("S", 0)])); // Trp->Arg
    composition.insert("[UNIMOD:676]".to_string(), HashMap::from([("H", -7), ("C", -9), ("N", -1), ("O", 0), ("S", 0)])); // Trp->Gly
    composition.insert("[UNIMOD:677]".to_string(), HashMap::from([("H", 1), ("C", -5), ("N", -1)])); // Trp->Xle
    composition.insert("[UNIMOD:678]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", -1), ("S", 0)])); // Tyr->Phe
    composition.insert("[UNIMOD:679]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", 0), ("S", 0)])); // Tyr->Ser
    composition.insert("[UNIMOD:680]".to_string(), HashMap::from([("H", -3), ("C", -5), ("N", 1), ("O", 0), ("S", 0)])); // Tyr->Asn
    composition.insert("[UNIMOD:681]".to_string(), HashMap::from([("H", -2), ("C", -3), ("N", 2), ("O", -1), ("S", 0)])); // Tyr->His
    composition.insert("[UNIMOD:682]".to_string(), HashMap::from([("H", -4), ("C", -5), ("N", 0), ("O", 1), ("S", 0)])); // Tyr->Asp
    composition.insert("[UNIMOD:683]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", -1), ("S", 1)])); // Tyr->Cys
    composition.insert("[UNIMOD:684]".to_string(), HashMap::from([("H", 12), ("C", 11), ("N", 1), ("O", 1), ("Br", 1)])); // BDMAPP
    composition.insert("[UNIMOD:685]".to_string(), HashMap::from([("H", 31), ("C", 18), ("N", 1), ("O", 4)])); // NA-LNO2
    composition.insert("[UNIMOD:686]".to_string(), HashMap::from([("H", 33), ("C", 18), ("N", 1), ("O", 4)])); // NA-OA-NO2
    composition.insert("[UNIMOD:687]".to_string(), HashMap::from([("H", -1), ("2H", 4), ("C", 6), ("N", 1), ("O", 1)])); // ICPL:2H(4)
    composition.insert("[UNIMOD:695]".to_string(), HashMap::from([("C", -6), ("13C", 6), ("N", -1), ("15N", 1)])); // Label:13C(6)15N(1)
    composition.insert("[UNIMOD:696]".to_string(), HashMap::from([("H", -9), ("2H", 9), ("C", -6), ("13C", 6), ("N", -2), ("15N", 2)])); // Label:2H(9)13C(6)15N(2)
    composition.insert("[UNIMOD:697]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 1), ("O", 1)])); // NIC
    composition.insert("[UNIMOD:698]".to_string(), HashMap::from([("H", 1), ("2H", 3), ("C", 6), ("N", 1), ("O", 1)])); // dNIC
    composition.insert("[UNIMOD:720]".to_string(), HashMap::from([("H", 14), ("C", 9), ("O", 1)])); // HNE-Delta:H(2)O
    composition.insert("[UNIMOD:721]".to_string(), HashMap::from([("H", 14), ("C", 9), ("O", 2)])); // 4-ONE
    composition.insert("[UNIMOD:723]".to_string(), HashMap::from([("H", 5), ("C", 2), ("O", 3), ("P", 1)])); // O-Dimethylphosphate
    composition.insert("[UNIMOD:724]".to_string(), HashMap::from([("H", 3), ("C", 1), ("O", 3), ("P", 1)])); // O-Methylphosphate
    composition.insert("[UNIMOD:725]".to_string(), HashMap::from([("H", 9), ("C", 4), ("O", 3), ("P", 1)])); // Diethylphosphate
    composition.insert("[UNIMOD:726]".to_string(), HashMap::from([("H", 5), ("C", 2), ("O", 3), ("P", 1)])); // Ethylphosphate
    composition.insert("[UNIMOD:727]".to_string(), HashMap::from([("H", 15), ("C", 7), ("O", 2), ("P", 1)])); // O-pinacolylmethylphosphonate
    composition.insert("[UNIMOD:728]".to_string(), HashMap::from([("H", 3), ("C", 1), ("O", 2), ("P", 1)])); // Methylphosphonate
    composition.insert("[UNIMOD:729]".to_string(), HashMap::from([("H", 9), ("C", 4), ("O", 2), ("P", 1)])); // O-Isopropylmethylphosphonate
    composition.insert("[UNIMOD:730]".to_string(), HashMap::from([("H", 24), ("C", 7), ("13C", 7), ("N", 3), ("15N", 1), ("O", 3)])); // iTRAQ8plex
    composition.insert("[UNIMOD:731]".to_string(), HashMap::from([("H", 24), ("C", 8), ("13C", 6), ("N", 2), ("15N", 2), ("O", 3)])); // iTRAQ8plex:13C(6)15N(2)
    composition.insert("[UNIMOD:734]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 1)])); // Ethanolamine
    composition.insert("[UNIMOD:735]".to_string(), HashMap::from([("H", 8), ("C", 4), ("O", 1), ("S", 2)])); // BEMAD_ST
    composition.insert("[UNIMOD:736]".to_string(), HashMap::from([("H", 8), ("C", 4), ("O", 2), ("S", 1)])); // BEMAD_C
    composition.insert("[UNIMOD:737]".to_string(), HashMap::from([("H", 20), ("C", 8), ("13C", 4), ("N", 1), ("15N", 1), ("O", 2)])); // TMT6plex
    composition.insert("[UNIMOD:738]".to_string(), HashMap::from([("H", 20), ("C", 11), ("13C", 1), ("N", 2), ("O", 2)])); // TMT2plex
    composition.insert("[UNIMOD:739]".to_string(), HashMap::from([("H", 20), ("C", 12), ("N", 2), ("O", 2)])); // TMT
    composition.insert("[UNIMOD:740]".to_string(), HashMap::from([("H", 50), ("C", 23), ("13C", 12), ("N", 8), ("15N", 6), ("O", 18)])); // ExacTagThiol
    composition.insert("[UNIMOD:741]".to_string(), HashMap::from([("H", 52), ("C", 25), ("13C", 12), ("N", 8), ("15N", 6), ("O", 19), ("S", 1)])); // ExacTagAmine
    composition.insert("[UNIMOD:743]".to_string(), HashMap::from([("H", 12), ("C", 9), ("O", 1)])); // 4-ONE+Delta:H(-2)O(-1)
    composition.insert("[UNIMOD:744]".to_string(), HashMap::from([("H", 9), ("C", 10), ("N", 3), ("O", 3), ("S", 1)])); // NO_SMX_SEMD
    composition.insert("[UNIMOD:746]".to_string(), HashMap::from([("H", 9), ("C", 10), ("N", 3), ("O", 4), ("S", 1)])); // NO_SMX_SIMD
    composition.insert("[UNIMOD:747]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 3)])); // Malonyl
    composition.insert("[UNIMOD:748]".to_string(), HashMap::from([("H", 4), ("C", 7), ("O", 4), ("S", 1)])); // 3sulfo
    composition.insert("[UNIMOD:750]".to_string(), HashMap::from([("H", -3), ("F", 3)])); // trifluoro
    composition.insert("[UNIMOD:751]".to_string(), HashMap::from([("H", 1), ("C", 6), ("N", 3), ("O", 6)])); // TNBS
    composition.insert("[UNIMOD:762]".to_string(), HashMap::from([("H", 7), ("C", 9), ("N", 1), ("O", 1), ("Cl", 2)])); // IDEnT
    composition.insert("[UNIMOD:763]".to_string(), HashMap::from([("H", 2), ("2H", 6), ("C", 4), ("O", 1), ("S", 2)])); // BEMAD_ST:2H(6)
    composition.insert("[UNIMOD:764]".to_string(), HashMap::from([("H", 2), ("2H", 6), ("C", 4), ("O", 2), ("S", 1)])); // BEMAD_C:2H(6)
    composition.insert("[UNIMOD:765]".to_string(), HashMap::from([("H", -9), ("C", -5), ("N", -1), ("O", -1), ("S", -1)])); // Met-loss
    composition.insert("[UNIMOD:766]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -1), ("S", -1)])); // Met-loss+Acetyl
    composition.insert("[UNIMOD:767]".to_string(), HashMap::from([("H", 8), ("C", 11), ("O", 2)])); // Menadione-HQ
    composition.insert("[UNIMOD:768]".to_string(), HashMap::from([("H", 1), ("2H", 3), ("C", 3), ("O", 1)])); // Methyl+Acetyl:2H(3)
    composition.insert("[UNIMOD:771]".to_string(), HashMap::from([("H", 16), ("C", 16), ("O", 2)])); // lapachenole
    composition.insert("[UNIMOD:772]".to_string(), HashMap::from([("C", -5), ("13C", 5)])); // Label:13C(5)
    composition.insert("[UNIMOD:773]".to_string(), HashMap::from([("H", 3), ("C", 4), ("N", 1), ("O", 2)])); // maleimide
    composition.insert("[UNIMOD:774]".to_string(), HashMap::from([("H", 38), ("C", 29), ("N", 8), ("O", 6), ("S", 1)])); // Biotin-phenacyl
    composition.insert("[UNIMOD:775]".to_string(), HashMap::from([("H", 2), ("13C", 2), ("O", 2)])); // Carboxymethyl:13C(2)
    composition.insert("[UNIMOD:776]".to_string(), HashMap::from([("H", 2), ("2H", 5), ("C", 6), ("N", 1), ("O", 2)])); // NEM:2H(5)
    composition.insert("[UNIMOD:792]".to_string(), HashMap::from([("H", 1), ("2H", 4), ("C", 2), ("N", 1), ("O", -1), ("S", 1)])); // AEC-MAEC:2H(4)
    composition.insert("[UNIMOD:793]".to_string(), HashMap::from([("H", 23), ("C", 14), ("N", 1), ("O", 10)])); // Hex(1)HexNAc(1)
    composition.insert("[UNIMOD:799]".to_string(), HashMap::from([("H", 6), ("C", -2), ("13C", 6), ("N", 2), ("O", 2)])); // Label:13C(6)+GG
    composition.insert("[UNIMOD:800]".to_string(), HashMap::from([("H", 25), ("C", 15), ("N", 3), ("O", 2), ("S", 1)])); // Biotin:Thermo-21345
    composition.insert("[UNIMOD:801]".to_string(), HashMap::from([("H", 10), ("C", 5)])); // Pentylamine
    composition.insert("[UNIMOD:811]".to_string(), HashMap::from([("H", 37), ("C", 21), ("N", 5), ("O", 6), ("S", 1)])); // Biotin:Thermo-21360
    composition.insert("[UNIMOD:821]".to_string(), HashMap::from([("H", 38), ("C", 37), ("N", 4), ("O", 7), ("S", 1)])); // Cy3b-maleimide
    composition.insert("[UNIMOD:822]".to_string(), HashMap::from([("H", -2), ("C", -2), ("O", -2)])); // Gly-loss+Amide
    composition.insert("[UNIMOD:824]".to_string(), HashMap::from([("H", 8), ("C", 10), ("N", 2), ("O", 4)])); // Xlink:BMOE
    composition.insert("[UNIMOD:825]".to_string(), HashMap::from([("C", 6), ("N", 2), ("O", 4)])); // Xlink:DFDNB
    composition.insert("[UNIMOD:827]".to_string(), HashMap::from([("H", 33), ("C", 29), ("O", 10), ("P", 1)])); // TMPP-Ac
    composition.insert("[UNIMOD:830]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 2)])); // Dihydroxyimidazolidine
    composition.insert("[UNIMOD:834]".to_string(), HashMap::from([("H", -2), ("2H", 4), ("C", 2), ("O", 1)])); // Label:2H(4)+Acetyl
    composition.insert("[UNIMOD:835]".to_string(), HashMap::from([("H", 2), ("C", -4), ("13C", 6), ("O", 1)])); // Label:13C(6)+Acetyl
    composition.insert("[UNIMOD:836]".to_string(), HashMap::from([("H", 2), ("C", -4), ("13C", 6), ("N", -2), ("15N", 2), ("O", 1)])); // Label:13C(6)15N(2)+Acetyl
    composition.insert("[UNIMOD:837]".to_string(), HashMap::from([("H", -1), ("C", 3), ("N", 1), ("O", 2)])); // Arg->Npo
    composition.insert("[UNIMOD:846]".to_string(), HashMap::from([("H", 32), ("C", 20), ("N", 6), ("O", 8)])); // EQIGG
    composition.insert("[UNIMOD:848]".to_string(), HashMap::from([("H", 10), ("C", 16), ("O", 4)])); // Arg2PG
    composition.insert("[UNIMOD:849]".to_string(), HashMap::from([("H", 10), ("C", 10), ("N", 5), ("O", 7), ("P", 1)])); // cGMP
    composition.insert("[UNIMOD:851]".to_string(), HashMap::from([("H", 4), ("C", 5), ("N", 5), ("O", 1)])); // cGMP+RMP-loss
    composition.insert("[UNIMOD:853]".to_string(), HashMap::from([("H", 2), ("2H", 4), ("C", 4), ("N", 2), ("O", 2)])); // Label:2H(4)+GG
    composition.insert("[UNIMOD:859]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 1)])); // MG-H1
    composition.insert("[UNIMOD:860]".to_string(), HashMap::from([("C", 2), ("O", 1)])); // G-H1
    composition.insert("[UNIMOD:861]".to_string(), HashMap::from([("H", 53), ("C", 37), ("N", 6), ("O", 6), ("F", 2), ("S", 1), ("B", 1)])); // ZGB
    composition.insert("[UNIMOD:862]".to_string(), HashMap::from([("H", -3), ("2H", 3), ("C", -1), ("13C", 1)])); // Label:13C(1)2H(3)
    composition.insert("[UNIMOD:864]".to_string(), HashMap::from([("H", 6), ("C", -2), ("13C", 6), ("15N", 2), ("O", 2)])); // Label:13C(6)15N(2)+GG
    composition.insert("[UNIMOD:866]".to_string(), HashMap::from([("H", -1), ("2H", 4), ("13C", 6), ("N", 1), ("O", 1)])); // ICPL:13C(6)2H(4)
    composition.insert("[UNIMOD:876]".to_string(), HashMap::from([("H", 36), ("C", 23), ("N", 8), ("O", 11)])); // QEQTGG
    composition.insert("[UNIMOD:877]".to_string(), HashMap::from([("H", 37), ("C", 23), ("N", 9), ("O", 10)])); // QQQTGG
    composition.insert("[UNIMOD:884]".to_string(), HashMap::from([("H", 45), ("C", 34), ("N", 7), ("O", 7), ("S", 1)])); // Biotin:Thermo-21325
    composition.insert("[UNIMOD:885]".to_string(), HashMap::from([("H", -3), ("2H", 3), ("C", -1), ("13C", 1), ("O", 1)])); // Label:13C(1)2H(3)+Oxidation
    composition.insert("[UNIMOD:886]".to_string(), HashMap::from([("H", 4), ("C", 6), ("O", 2)])); // HydroxymethylOP
    composition.insert("[UNIMOD:887]".to_string(), HashMap::from([("H", 21), ("C", 20), ("N", 3), ("O", 5)])); // MDCC
    composition.insert("[UNIMOD:888]".to_string(), HashMap::from([("H", 12), ("C", 7), ("N", 2), ("O", 1)])); // mTRAQ
    composition.insert("[UNIMOD:889]".to_string(), HashMap::from([("H", 12), ("C", 4), ("13C", 3), ("N", 1), ("15N", 1), ("O", 1)])); // mTRAQ:13C(3)15N(1)
    composition.insert("[UNIMOD:890]".to_string(), HashMap::from([("H", 48), ("C", 39), ("N", 4), ("O", 15), ("S", 4)])); // DyLight-maleimide
    composition.insert("[UNIMOD:891]".to_string(), HashMap::from([("H", 58), ("C", 32), ("N", 2), ("O", 15)])); // Methyl-PEO12-Maleimide
    composition.insert("[UNIMOD:893]".to_string(), HashMap::from([("H", 11), ("C", 6), ("N", 1), ("O", 3), ("S", 2)])); // CarbamidomethylDTT
    composition.insert("[UNIMOD:894]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 4), ("S", 2)])); // CarboxymethylDTT
    composition.insert("[UNIMOD:895]".to_string(), HashMap::from([("H", 42), ("C", 26), ("N", 8), ("O", 7)])); // Biotin-PEG-PRA
    composition.insert("[UNIMOD:896]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", 3), ("S", -1)])); // Met->Aha
    composition.insert("[UNIMOD:897]".to_string(), HashMap::from([("N", -4), ("15N", 4)])); // Label:15N(4)
    composition.insert("[UNIMOD:898]".to_string(), HashMap::from([("H", 2), ("O", 6), ("P", 2)])); // pyrophospho
    composition.insert("[UNIMOD:899]".to_string(), HashMap::from([("H", -2), ("C", 1), ("S", -1)])); // Met->Hpg
    composition.insert("[UNIMOD:901]".to_string(), HashMap::from([("H", 24), ("C", 17), ("O", 9)])); // 4AcAllylGal
    composition.insert("[UNIMOD:902]".to_string(), HashMap::from([("H", 5), ("C", 2), ("As", 1)])); // DimethylArsino
    composition.insert("[UNIMOD:903]".to_string(), HashMap::from([("H", -4), ("C", -1), ("O", 1), ("S", 1)])); // Lys->CamCys
    composition.insert("[UNIMOD:904]".to_string(), HashMap::from([("H", -1), ("C", -4), ("N", 1), ("O", 1), ("S", 1)])); // Phe->CamCys
    composition.insert("[UNIMOD:905]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", 1), ("S", 1)])); // Leu->MetOx
    composition.insert("[UNIMOD:906]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -1), ("O", 1), ("S", 1)])); // Lys->MetOx
    composition.insert("[UNIMOD:907]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 6)])); // Galactosyl
    composition.insert("[UNIMOD:908]".to_string(), HashMap::from([("H", 27), ("C", 17), ("N", 3), ("O", 3)])); // Xlink:SMCC[321]
    composition.insert("[UNIMOD:910]".to_string(), HashMap::from([("H", 16), ("C", 10), ("N", 2), ("O", 4)])); // Bacillosamine
    composition.insert("[UNIMOD:911]".to_string(), HashMap::from([("H", 14), ("C", 9), ("N", 1), ("O", 1), ("S", 1)])); // MTSL
    composition.insert("[UNIMOD:912]".to_string(), HashMap::from([("H", 45), ("C", 25), ("N", 5), ("O", 4), ("S", 1)])); // HNE-BAHAH
    composition.insert("[UNIMOD:914]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 3)])); // Methylmalonylation
    composition.insert("[UNIMOD:923]".to_string(), HashMap::from([("H", 6), ("13C", 4), ("15N", 2), ("O", 2)])); // Label:13C(4)15N(2)+GG
    composition.insert("[UNIMOD:926]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 1), ("O", -1)])); // ethylamino
    composition.insert("[UNIMOD:928]".to_string(), HashMap::from([("H", 4), ("C", 2), ("S", 1)])); // MercaptoEthanol
    composition.insert("[UNIMOD:931]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", -1), ("O", 1)])); // Ethyl+Deamidated
    composition.insert("[UNIMOD:932]".to_string(), HashMap::from([("H", 55), ("C", 37), ("N", 11), ("O", 12)])); // VFQQQTGG
    composition.insert("[UNIMOD:933]".to_string(), HashMap::from([("H", 81), ("C", 53), ("N", 13), ("O", 19)])); // VIEVYQEQTGG
    composition.insert("[UNIMOD:934]".to_string(), HashMap::from([("H", 30), ("C", 19), ("N", 6), ("O", 10)])); // AMTzHexNAc2
    composition.insert("[UNIMOD:935]".to_string(), HashMap::from([("H", 32), ("C", 27), ("N", 5), ("O", 3)])); // Atto495Maleimide
    composition.insert("[UNIMOD:936]".to_string(), HashMap::from([("H", -1), ("Cl", 1)])); // Chlorination
    composition.insert("[UNIMOD:937]".to_string(), HashMap::from([("H", -2), ("Cl", 2)])); // dichlorination
    composition.insert("[UNIMOD:938]".to_string(), HashMap::from([("H", 52), ("C", 35), ("N", 10), ("O", 9), ("S", 2)])); // AROD
    composition.insert("[UNIMOD:939]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 1), ("S", -1)])); // Cys->methylaminoAla
    composition.insert("[UNIMOD:940]".to_string(), HashMap::from([("H", 5), ("C", 2), ("N", 1), ("S", -1)])); // Cys->ethylaminoAla
    composition.insert("[UNIMOD:941]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 2), ("O", 4), ("S", 1)])); // DNPS
    composition.insert("[UNIMOD:942]".to_string(), HashMap::from([("H", 26), ("C", 22), ("N", 4), ("O", 5), ("S", 1)])); // SulfoGMBS
    composition.insert("[UNIMOD:943]".to_string(), HashMap::from([("H", 21), ("C", 13), ("N", 3), ("O", 3)])); // DimethylamineGMBS
    composition.insert("[UNIMOD:944]".to_string(), HashMap::from([("H", -9), ("2H", 9), ("N", -2), ("15N", 2)])); // Label:15N(2)2H(9)
    composition.insert("[UNIMOD:946]".to_string(), HashMap::from([("H", 26), ("C", 20), ("O", 3)])); // LG-anhydrolactam
    composition.insert("[UNIMOD:947]".to_string(), HashMap::from([("H", 28), ("C", 20), ("O", 3)])); // LG-pyrrole
    composition.insert("[UNIMOD:948]".to_string(), HashMap::from([("H", 26), ("C", 20), ("O", 2)])); // LG-anhyropyrrole
    composition.insert("[UNIMOD:949]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 4)])); // 3-deoxyglucosone
    composition.insert("[UNIMOD:950]".to_string(), HashMap::from([("H", -1), ("Li", 1)])); // Cation:Li
    composition.insert("[UNIMOD:951]".to_string(), HashMap::from([("H", -2), ("Ca", 1)])); // Cation:Ca[II]
    composition.insert("[UNIMOD:952]".to_string(), HashMap::from([("H", -2), ("Fe", 1)])); // Cation:Fe[II]
    composition.insert("[UNIMOD:953]".to_string(), HashMap::from([("H", -2), ("Ni", 1)])); // Cation:Ni[II]
    composition.insert("[UNIMOD:954]".to_string(), HashMap::from([("H", -2), ("Zn", 1)])); // Cation:Zn[II]
    composition.insert("[UNIMOD:955]".to_string(), HashMap::from([("H", -1), ("Ag", 1)])); // Cation:Ag
    composition.insert("[UNIMOD:956]".to_string(), HashMap::from([("H", -2), ("Mg", 1)])); // Cation:Mg[II]
    composition.insert("[UNIMOD:957]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 4)])); // 2-succinyl
    composition.insert("[UNIMOD:958]".to_string(), HashMap::from([("H", 3), ("C", 3), ("N", 1), ("O", -1)])); // Propargylamine
    composition.insert("[UNIMOD:959]".to_string(), HashMap::from([("H", 4), ("C", 3), ("N", 1), ("O", 2), ("P", 1)])); // Phosphopropargyl
    composition.insert("[UNIMOD:960]".to_string(), HashMap::from([("H", 137), ("C", 90), ("N", 21), ("O", 37), ("S", 1)])); // SUMO2135
    composition.insert("[UNIMOD:961]".to_string(), HashMap::from([("H", 224), ("C", 150), ("N", 38), ("O", 60), ("S", 1)])); // SUMO3549
    composition.insert("[UNIMOD:967]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 2), ("S", 1)])); // thioacylPA
    composition.insert("[UNIMOD:971]".to_string(), HashMap::from([("H", 59), ("C", 37), ("N", 7), ("O", 23)])); // maleimide3
    composition.insert("[UNIMOD:972]".to_string(), HashMap::from([("H", 79), ("C", 49), ("N", 7), ("O", 33)])); // maleimide5
    composition.insert("[UNIMOD:973]".to_string(), HashMap::from([("H", 27), ("C", 22), ("N", 7), ("O", 4)])); // Puromycin
    composition.insert("[UNIMOD:977]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1)])); // Carbofuran
    composition.insert("[UNIMOD:978]".to_string(), HashMap::from([("H", 7), ("C", 8), ("N", 1), ("S", 1)])); // BITC
    composition.insert("[UNIMOD:979]".to_string(), HashMap::from([("H", 9), ("C", 9), ("N", 1), ("S", 1)])); // PEITC
    composition.insert("[UNIMOD:981]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 5)])); // glucosone
    composition.insert("[UNIMOD:984]".to_string(), HashMap::from([("H", 25), ("C", 14), ("N", 3), ("O", 2), ("S", 1)])); // cysTMT
    composition.insert("[UNIMOD:985]".to_string(), HashMap::from([("H", 25), ("C", 10), ("13C", 4), ("N", 2), ("15N", 1), ("O", 2), ("S", 1)])); // cysTMT6plex
    composition.insert("[UNIMOD:986]".to_string(), HashMap::from([("H", 4), ("C", -4), ("13C", 6)])); // Label:13C(6)+Dimethyl
    composition.insert("[UNIMOD:987]".to_string(), HashMap::from([("H", 4), ("C", -4), ("13C", 6), ("N", -2), ("15N", 2)])); // Label:13C(6)15N(2)+Dimethyl
    composition.insert("[UNIMOD:989]".to_string(), HashMap::from([("H", 3), ("N", 1)])); // Ammonium
    composition.insert("[UNIMOD:991]".to_string(), HashMap::from([("H", -1), ("N", -1)])); // ISD_z+2_ion
    composition.insert("[UNIMOD:993]".to_string(), HashMap::from([("H", 27), ("C", 20), ("N", 5), ("O", 5), ("S", 1)])); // Biotin:Sigma-B1267
    composition.insert("[UNIMOD:994]".to_string(), HashMap::from([("N", -1), ("15N", 1)])); // Label:15N(1)
    composition.insert("[UNIMOD:995]".to_string(), HashMap::from([("N", -2), ("15N", 2)])); // Label:15N(2)
    composition.insert("[UNIMOD:996]".to_string(), HashMap::from([("N", -3), ("15N", 3)])); // Label:15N(3)
    composition.insert("[UNIMOD:997]".to_string(), HashMap::from([("H", 1), ("N", 1), ("O", 3), ("S", 1)])); // sulfo+amino
    composition.insert("[UNIMOD:1000]".to_string(), HashMap::from([("H", 5), ("C", 4), ("N", 5), ("O", 1), ("S", -1)])); // AHA-Alkyne
    composition.insert("[UNIMOD:1001]".to_string(), HashMap::from([("H", 37), ("C", 26), ("N", 11), ("O", 14), ("S", -1)])); // AHA-Alkyne-KDDDD
    composition.insert("[UNIMOD:1002]".to_string(), HashMap::from([("H", 16), ("C", 22), ("O", 11)])); // EGCG1
    composition.insert("[UNIMOD:1003]".to_string(), HashMap::from([("H", 11), ("C", 15), ("O", 6)])); // EGCG2
    composition.insert("[UNIMOD:1004]".to_string(), HashMap::from([("H", 2), ("C", -5), ("13C", 6), ("N", -4), ("15N", 4)])); // Label:13C(6)15N(4)+Methyl
    composition.insert("[UNIMOD:1005]".to_string(), HashMap::from([("H", 4), ("C", -4), ("13C", 6), ("N", -4), ("15N", 4)])); // Label:13C(6)15N(4)+Dimethyl
    composition.insert("[UNIMOD:1006]".to_string(), HashMap::from([("H", -1), ("2H", 3), ("C", -6), ("13C", 7), ("N", -4), ("15N", 4)])); // Label:13C(6)15N(4)+Methyl:2H(3)13C(1)
    composition.insert("[UNIMOD:1007]".to_string(), HashMap::from([("H", -2), ("2H", 6), ("C", -6), ("13C", 8), ("N", -4), ("15N", 4)])); // Label:13C(6)15N(4)+Dimethyl:2H(6)13C(2)
    composition.insert("[UNIMOD:1008]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1), ("S", -1), ("Se", 1)])); // Cys->CamSec
    composition.insert("[UNIMOD:1009]".to_string(), HashMap::from([("C", 1)])); // Thiazolidine
    composition.insert("[UNIMOD:1010]".to_string(), HashMap::from([("H", 122), ("C", 89), ("N", 18), ("O", 31), ("S", 1)])); // DEDGFLYMVYASQETFG
    composition.insert("[UNIMOD:1012]".to_string(), HashMap::from([("H", 33), ("C", 23), ("N", 5), ("O", 7), ("S", 1)])); // Biotin:Invitrogen-M1602
    composition.insert("[UNIMOD:1014]".to_string(), HashMap::from([("H", 5), ("C", 3), ("N", 1), ("O", 2)])); // glycidamide
    composition.insert("[UNIMOD:1015]".to_string(), HashMap::from([("H", 27), ("C", 16), ("N", 3), ("O", 3)])); // Ahx2+Hsl
    composition.insert("[UNIMOD:1017]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 1)])); // DMPO
    composition.insert("[UNIMOD:1018]".to_string(), HashMap::from([("H", 10), ("C", 8), ("O", 2)])); // ICDID
    composition.insert("[UNIMOD:1019]".to_string(), HashMap::from([("H", 4), ("2H", 6), ("C", 8), ("O", 2)])); // ICDID:2H(6)
    composition.insert("[UNIMOD:1020]".to_string(), HashMap::from([("H", 12), ("C", 8), ("O", 3)])); // Xlink:DSS[156]
    composition.insert("[UNIMOD:1021]".to_string(), HashMap::from([("H", 12), ("C", 10), ("O", 7)])); // Xlink:EGS[244]
    composition.insert("[UNIMOD:1022]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 5)])); // Xlink:DST[132]
    composition.insert("[UNIMOD:1023]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 3), ("S", 2)])); // Xlink:DTSSP[192]
    composition.insert("[UNIMOD:1024]".to_string(), HashMap::from([("H", 15), ("C", 12), ("N", 1), ("O", 4)])); // Xlink:SMCC[237]
    composition.insert("[UNIMOD:1027]".to_string(), HashMap::from([("H", 12), ("C", 7), ("N", 2), ("O", 1)])); // Xlink:DMP[140]
    composition.insert("[UNIMOD:1028]".to_string(), HashMap::from([("H", 5), ("C", 4), ("N", 1), ("O", 3)])); // Xlink:EGS[115]
    composition.insert("[UNIMOD:1031]".to_string(), HashMap::from([("H", 16), ("C", 10), ("N", 2), ("O", 2)])); // Biotin:Thermo-88310
    composition.insert("[UNIMOD:1032]".to_string(), HashMap::from([("H", 5), ("C", 7), ("N", 1), ("O", 2)])); // 2-nitrobenzyl
    composition.insert("[UNIMOD:1033]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 1), ("O", 2), ("S", -1), ("Se", 1)])); // Cys->SecNEM
    composition.insert("[UNIMOD:1034]".to_string(), HashMap::from([("H", 2), ("2H", 5), ("C", 6), ("N", 1), ("O", 2), ("S", -1), ("Se", 1)])); // Cys->SecNEM:2H(5)
    composition.insert("[UNIMOD:1035]".to_string(), HashMap::from([("H", 6), ("C", 9), ("N", 2), ("S", 1)])); // Thiadiazole
    composition.insert("[UNIMOD:1036]".to_string(), HashMap::from([("H", 38), ("C", 28), ("O", 6)])); // Withaferin
    composition.insert("[UNIMOD:1037]".to_string(), HashMap::from([("H", 42), ("C", 22), ("N", 3), ("O", 4), ("P", 1)])); // Biotin:Thermo-88317
    composition.insert("[UNIMOD:1038]".to_string(), HashMap::from([("H", 46), ("C", 37), ("N", 3), ("O", 6), ("P", 1)])); // TAMRA-FP
    composition.insert("[UNIMOD:1039]".to_string(), HashMap::from([("H", 37), ("C", 23), ("N", 5), ("O", 8), ("S", 1)])); // Biotin:Thermo-21901+H2O
    composition.insert("[UNIMOD:1041]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 1)])); // Deoxyhypusine
    composition.insert("[UNIMOD:1042]".to_string(), HashMap::from([("H", 11), ("C", 6), ("N", 1), ("O", 1)])); // Acetyldeoxyhypusine
    composition.insert("[UNIMOD:1043]".to_string(), HashMap::from([("H", 11), ("C", 6), ("N", 1), ("O", 2)])); // Acetylhypusine
    composition.insert("[UNIMOD:1044]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 0), ("S", 1)])); // Ala->Cys
    composition.insert("[UNIMOD:1045]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", 0), ("S", 0)])); // Ala->Phe
    composition.insert("[UNIMOD:1046]".to_string(), HashMap::from([("H", 2), ("C", 3), ("N", 2), ("O", 0), ("S", 0)])); // Ala->His
    composition.insert("[UNIMOD:1047]".to_string(), HashMap::from([("H", 6), ("C", 3)])); // Ala->Xle
    composition.insert("[UNIMOD:1048]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 1), ("O", 0), ("S", 0)])); // Ala->Lys
    composition.insert("[UNIMOD:1049]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", 0), ("S", 1)])); // Ala->Met
    composition.insert("[UNIMOD:1050]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", 1), ("O", 1), ("S", 0)])); // Ala->Asn
    composition.insert("[UNIMOD:1051]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1), ("S", 0)])); // Ala->Gln
    composition.insert("[UNIMOD:1052]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 3), ("O", 0), ("S", 0)])); // Ala->Arg
    composition.insert("[UNIMOD:1053]".to_string(), HashMap::from([("H", 5), ("C", 8), ("N", 1), ("O", 0), ("S", 0)])); // Ala->Trp
    composition.insert("[UNIMOD:1054]".to_string(), HashMap::from([("H", 4), ("C", 6), ("N", 0), ("O", 1), ("S", 0)])); // Ala->Tyr
    composition.insert("[UNIMOD:1055]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 0), ("S", -1)])); // Cys->Ala
    composition.insert("[UNIMOD:1056]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 0), ("O", 2), ("S", -1)])); // Cys->Asp
    composition.insert("[UNIMOD:1057]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 2), ("S", -1)])); // Cys->Glu
    composition.insert("[UNIMOD:1058]".to_string(), HashMap::from([("H", 2), ("C", 3), ("N", 2), ("O", 0), ("S", -1)])); // Cys->His
    composition.insert("[UNIMOD:1059]".to_string(), HashMap::from([("H", 6), ("C", 3), ("S", -1)])); // Cys->Xle
    composition.insert("[UNIMOD:1060]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 1), ("O", 0), ("S", -1)])); // Cys->Lys
    composition.insert("[UNIMOD:1061]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", 0), ("S", 0)])); // Cys->Met
    composition.insert("[UNIMOD:1062]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", 1), ("O", 1), ("S", -1)])); // Cys->Asn
    composition.insert("[UNIMOD:1063]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 0), ("S", -1)])); // Cys->Pro
    composition.insert("[UNIMOD:1064]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1), ("S", -1)])); // Cys->Gln
    composition.insert("[UNIMOD:1065]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 1), ("S", -1)])); // Cys->Thr
    composition.insert("[UNIMOD:1066]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", 0), ("S", -1)])); // Cys->Val
    composition.insert("[UNIMOD:1067]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", 0), ("O", -2), ("S", 1)])); // Asp->Cys
    composition.insert("[UNIMOD:1068]".to_string(), HashMap::from([("H", 4), ("C", 5), ("N", 0), ("O", -2), ("S", 0)])); // Asp->Phe
    composition.insert("[UNIMOD:1069]".to_string(), HashMap::from([("H", 6), ("C", 2), ("O", -2)])); // Asp->Xle
    composition.insert("[UNIMOD:1070]".to_string(), HashMap::from([("H", 7), ("C", 2), ("N", 1), ("O", -2), ("S", 0)])); // Asp->Lys
    composition.insert("[UNIMOD:1071]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 0), ("O", -2), ("S", 1)])); // Asp->Met
    composition.insert("[UNIMOD:1072]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", -2), ("S", 0)])); // Asp->Pro
    composition.insert("[UNIMOD:1073]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 1), ("O", -1), ("S", 0)])); // Asp->Gln
    composition.insert("[UNIMOD:1074]".to_string(), HashMap::from([("H", 7), ("C", 2), ("N", 3), ("O", -2), ("S", 0)])); // Asp->Arg
    composition.insert("[UNIMOD:1075]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", 0), ("O", -1), ("S", 0)])); // Asp->Ser
    composition.insert("[UNIMOD:1076]".to_string(), HashMap::from([("H", 2), ("C", 0), ("N", 0), ("O", -1), ("S", 0)])); // Asp->Thr
    composition.insert("[UNIMOD:1077]".to_string(), HashMap::from([("H", 5), ("C", 7), ("N", 1), ("O", -2), ("S", 0)])); // Asp->Trp
    composition.insert("[UNIMOD:1078]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", -2), ("S", 1)])); // Glu->Cys
    composition.insert("[UNIMOD:1079]".to_string(), HashMap::from([("H", 2), ("C", 4), ("N", 0), ("O", -2), ("S", 0)])); // Glu->Phe
    composition.insert("[UNIMOD:1080]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 2), ("O", -2), ("S", 0)])); // Glu->His
    composition.insert("[UNIMOD:1081]".to_string(), HashMap::from([("H", 4), ("C", 1), ("O", -2)])); // Glu->Xle
    composition.insert("[UNIMOD:1082]".to_string(), HashMap::from([("H", 2), ("C", 0), ("N", 0), ("O", -2), ("S", 1)])); // Glu->Met
    composition.insert("[UNIMOD:1083]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", 1), ("O", -1), ("S", 0)])); // Glu->Asn
    composition.insert("[UNIMOD:1084]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", -2), ("S", 0)])); // Glu->Pro
    composition.insert("[UNIMOD:1085]".to_string(), HashMap::from([("H", 5), ("C", 1), ("N", 3), ("O", -2), ("S", 0)])); // Glu->Arg
    composition.insert("[UNIMOD:1086]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", -1), ("S", 0)])); // Glu->Ser
    composition.insert("[UNIMOD:1087]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", 0), ("O", -1), ("S", 0)])); // Glu->Thr
    composition.insert("[UNIMOD:1088]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 1), ("O", -2), ("S", 0)])); // Glu->Trp
    composition.insert("[UNIMOD:1089]".to_string(), HashMap::from([("H", 2), ("C", 4), ("N", 0), ("O", -1), ("S", 0)])); // Glu->Tyr
    composition.insert("[UNIMOD:1090]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", 0), ("S", 0)])); // Phe->Ala
    composition.insert("[UNIMOD:1091]".to_string(), HashMap::from([("H", -4), ("C", -5), ("N", 0), ("O", 2), ("S", 0)])); // Phe->Asp
    composition.insert("[UNIMOD:1092]".to_string(), HashMap::from([("H", -2), ("C", -4), ("N", 0), ("O", 2), ("S", 0)])); // Phe->Glu
    composition.insert("[UNIMOD:1093]".to_string(), HashMap::from([("H", -6), ("C", -7), ("N", 0), ("O", 0), ("S", 0)])); // Phe->Gly
    composition.insert("[UNIMOD:1094]".to_string(), HashMap::from([("H", -2), ("C", -3), ("N", 2), ("O", 0), ("S", 0)])); // Phe->His
    composition.insert("[UNIMOD:1095]".to_string(), HashMap::from([("H", 3), ("C", -3), ("N", 1), ("O", 0), ("S", 0)])); // Phe->Lys
    composition.insert("[UNIMOD:1096]".to_string(), HashMap::from([("H", 0), ("C", -4), ("N", 0), ("O", 0), ("S", 1)])); // Phe->Met
    composition.insert("[UNIMOD:1097]".to_string(), HashMap::from([("H", -3), ("C", -5), ("N", 1), ("O", 1), ("S", 0)])); // Phe->Asn
    composition.insert("[UNIMOD:1098]".to_string(), HashMap::from([("H", -2), ("C", -4), ("N", 0), ("O", 0), ("S", 0)])); // Phe->Pro
    composition.insert("[UNIMOD:1099]".to_string(), HashMap::from([("H", -1), ("C", -4), ("N", 1), ("O", 1), ("S", 0)])); // Phe->Gln
    composition.insert("[UNIMOD:1100]".to_string(), HashMap::from([("H", 3), ("C", -3), ("N", 3), ("O", 0), ("S", 0)])); // Phe->Arg
    composition.insert("[UNIMOD:1101]".to_string(), HashMap::from([("H", -2), ("C", -5), ("N", 0), ("O", 1), ("S", 0)])); // Phe->Thr
    composition.insert("[UNIMOD:1102]".to_string(), HashMap::from([("H", 1), ("C", 2), ("N", 1), ("O", 0), ("S", 0)])); // Phe->Trp
    composition.insert("[UNIMOD:1103]".to_string(), HashMap::from([("H", 6), ("C", 7), ("N", 0), ("O", 0), ("S", 0)])); // Gly->Phe
    composition.insert("[UNIMOD:1104]".to_string(), HashMap::from([("H", 4), ("C", 4), ("N", 2), ("O", 0), ("S", 0)])); // Gly->His
    composition.insert("[UNIMOD:1105]".to_string(), HashMap::from([("H", 8), ("C", 4)])); // Gly->Xle
    composition.insert("[UNIMOD:1106]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 1), ("O", 0), ("S", 0)])); // Gly->Lys
    composition.insert("[UNIMOD:1107]".to_string(), HashMap::from([("H", 6), ("C", 3), ("N", 0), ("O", 0), ("S", 1)])); // Gly->Met
    composition.insert("[UNIMOD:1108]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1), ("S", 0)])); // Gly->Asn
    composition.insert("[UNIMOD:1109]".to_string(), HashMap::from([("H", 4), ("C", 3), ("N", 0), ("O", 0), ("S", 0)])); // Gly->Pro
    composition.insert("[UNIMOD:1110]".to_string(), HashMap::from([("H", 5), ("C", 3), ("N", 1), ("O", 1), ("S", 0)])); // Gly->Gln
    composition.insert("[UNIMOD:1111]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", 1), ("S", 0)])); // Gly->Thr
    composition.insert("[UNIMOD:1112]".to_string(), HashMap::from([("H", 6), ("C", 7), ("N", 0), ("O", 1), ("S", 0)])); // Gly->Tyr
    composition.insert("[UNIMOD:1113]".to_string(), HashMap::from([("H", -2), ("C", -3), ("N", -2), ("O", 0), ("S", 0)])); // His->Ala
    composition.insert("[UNIMOD:1114]".to_string(), HashMap::from([("H", -2), ("C", -3), ("N", -2), ("O", 0), ("S", 1)])); // His->Cys
    composition.insert("[UNIMOD:1115]".to_string(), HashMap::from([("H", 0), ("C", -1), ("N", -2), ("O", 2), ("S", 0)])); // His->Glu
    composition.insert("[UNIMOD:1116]".to_string(), HashMap::from([("H", 2), ("C", 3), ("N", -2), ("O", 0), ("S", 0)])); // His->Phe
    composition.insert("[UNIMOD:1117]".to_string(), HashMap::from([("H", -4), ("C", -4), ("N", -2), ("O", 0), ("S", 0)])); // His->Gly
    composition.insert("[UNIMOD:1119]".to_string(), HashMap::from([("H", 5), ("C", 0), ("N", -1), ("O", 0), ("S", 0)])); // His->Lys
    composition.insert("[UNIMOD:1120]".to_string(), HashMap::from([("H", 2), ("C", -1), ("N", -2), ("O", 0), ("S", 1)])); // His->Met
    composition.insert("[UNIMOD:1121]".to_string(), HashMap::from([("H", -2), ("C", -3), ("N", -2), ("O", 1), ("S", 0)])); // His->Ser
    composition.insert("[UNIMOD:1122]".to_string(), HashMap::from([("H", 0), ("C", -2), ("N", -2), ("O", 1), ("S", 0)])); // His->Thr
    composition.insert("[UNIMOD:1123]".to_string(), HashMap::from([("H", 2), ("C", -1), ("N", -2), ("O", 0), ("S", 0)])); // His->Val
    composition.insert("[UNIMOD:1124]".to_string(), HashMap::from([("H", 3), ("C", 5), ("N", -1), ("O", 0), ("S", 0)])); // His->Trp
    composition.insert("[UNIMOD:1125]".to_string(), HashMap::from([("H", -6), ("C", -3), ("N", 0), ("O", 0), ("S", 0)])); // Xle->Ala
    composition.insert("[UNIMOD:1126]".to_string(), HashMap::from([("H", -6), ("C", -3), ("N", 0), ("O", 0), ("S", 1)])); // Xle->Cys
    composition.insert("[UNIMOD:1127]".to_string(), HashMap::from([("H", -6), ("C", -2), ("N", 0), ("O", 2), ("S", 0)])); // Xle->Asp
    composition.insert("[UNIMOD:1128]".to_string(), HashMap::from([("H", -4), ("C", -1), ("N", 0), ("O", 2), ("S", 0)])); // Xle->Glu
    composition.insert("[UNIMOD:1129]".to_string(), HashMap::from([("H", -8), ("C", -4), ("N", 0), ("O", 0), ("S", 0)])); // Xle->Gly
    composition.insert("[UNIMOD:1130]".to_string(), HashMap::from([("H", -2), ("C", 3), ("N", 0), ("O", 1), ("S", 0)])); // Xle->Tyr
    composition.insert("[UNIMOD:1131]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -1), ("O", 0), ("S", 0)])); // Lys->Ala
    composition.insert("[UNIMOD:1132]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -1), ("O", 0), ("S", 1)])); // Lys->Cys
    composition.insert("[UNIMOD:1133]".to_string(), HashMap::from([("H", -7), ("C", -2), ("N", -1), ("O", 2), ("S", 0)])); // Lys->Asp
    composition.insert("[UNIMOD:1134]".to_string(), HashMap::from([("H", -3), ("C", 3), ("N", -1), ("O", 0), ("S", 0)])); // Lys->Phe
    composition.insert("[UNIMOD:1135]".to_string(), HashMap::from([("H", -9), ("C", -4), ("N", -1), ("O", 0), ("S", 0)])); // Lys->Gly
    composition.insert("[UNIMOD:1136]".to_string(), HashMap::from([("H", -5), ("C", 0), ("N", 1), ("O", 0), ("S", 0)])); // Lys->His
    composition.insert("[UNIMOD:1137]".to_string(), HashMap::from([("H", -5), ("C", -1), ("N", -1), ("O", 0), ("S", 0)])); // Lys->Pro
    composition.insert("[UNIMOD:1138]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -1), ("O", 1), ("S", 0)])); // Lys->Ser
    composition.insert("[UNIMOD:1139]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -1), ("O", 0), ("S", 0)])); // Lys->Val
    composition.insert("[UNIMOD:1140]".to_string(), HashMap::from([("H", -2), ("C", 5), ("N", 0), ("O", 0), ("S", 0)])); // Lys->Trp
    composition.insert("[UNIMOD:1141]".to_string(), HashMap::from([("H", -3), ("C", 3), ("N", -1), ("O", 1), ("S", 0)])); // Lys->Tyr
    composition.insert("[UNIMOD:1142]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 0), ("S", -1)])); // Met->Ala
    composition.insert("[UNIMOD:1143]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 0), ("S", 0)])); // Met->Cys
    composition.insert("[UNIMOD:1144]".to_string(), HashMap::from([("H", -4), ("C", -1), ("N", 0), ("O", 2), ("S", -1)])); // Met->Asp
    composition.insert("[UNIMOD:1145]".to_string(), HashMap::from([("H", -2), ("C", 0), ("N", 0), ("O", 2), ("S", -1)])); // Met->Glu
    composition.insert("[UNIMOD:1146]".to_string(), HashMap::from([("H", 0), ("C", 4), ("N", 0), ("O", 0), ("S", -1)])); // Met->Phe
    composition.insert("[UNIMOD:1147]".to_string(), HashMap::from([("H", -6), ("C", -3), ("N", 0), ("O", 0), ("S", -1)])); // Met->Gly
    composition.insert("[UNIMOD:1148]".to_string(), HashMap::from([("H", -2), ("C", 1), ("N", 2), ("O", 0), ("S", -1)])); // Met->His
    composition.insert("[UNIMOD:1149]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", 1), ("O", 1), ("S", -1)])); // Met->Asn
    composition.insert("[UNIMOD:1150]".to_string(), HashMap::from([("H", -2), ("C", 0), ("N", 0), ("O", 0), ("S", -1)])); // Met->Pro
    composition.insert("[UNIMOD:1151]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", 1), ("O", 1), ("S", -1)])); // Met->Gln
    composition.insert("[UNIMOD:1152]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 1), ("S", -1)])); // Met->Ser
    composition.insert("[UNIMOD:1153]".to_string(), HashMap::from([("H", 1), ("C", 6), ("N", 1), ("O", 0), ("S", -1)])); // Met->Trp
    composition.insert("[UNIMOD:1154]".to_string(), HashMap::from([("H", 0), ("C", 4), ("N", 0), ("O", 1), ("S", -1)])); // Met->Tyr
    composition.insert("[UNIMOD:1155]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", -1), ("O", -1), ("S", 0)])); // Asn->Ala
    composition.insert("[UNIMOD:1156]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", -1), ("O", -1), ("S", 1)])); // Asn->Cys
    composition.insert("[UNIMOD:1157]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", -1), ("O", 1), ("S", 0)])); // Asn->Glu
    composition.insert("[UNIMOD:1158]".to_string(), HashMap::from([("H", 3), ("C", 5), ("N", -1), ("O", -1), ("S", 0)])); // Asn->Phe
    composition.insert("[UNIMOD:1159]".to_string(), HashMap::from([("H", -3), ("C", -2), ("N", -1), ("O", -1), ("S", 0)])); // Asn->Gly
    composition.insert("[UNIMOD:1160]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", -1), ("O", -1), ("S", 1)])); // Asn->Met
    composition.insert("[UNIMOD:1161]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", -1), ("O", -1), ("S", 0)])); // Asn->Pro
    composition.insert("[UNIMOD:1162]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", 0), ("S", 0)])); // Asn->Gln
    composition.insert("[UNIMOD:1163]".to_string(), HashMap::from([("H", 6), ("C", 2), ("N", 2), ("O", -1), ("S", 0)])); // Asn->Arg
    composition.insert("[UNIMOD:1164]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", -1), ("O", -1), ("S", 0)])); // Asn->Val
    composition.insert("[UNIMOD:1165]".to_string(), HashMap::from([("H", 4), ("C", 7), ("N", 0), ("O", -1), ("S", 0)])); // Asn->Trp
    composition.insert("[UNIMOD:1166]".to_string(), HashMap::from([("H", -2), ("C", -2), ("N", 0), ("O", 0), ("S", 1)])); // Pro->Cys
    composition.insert("[UNIMOD:1167]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 2), ("S", 0)])); // Pro->Asp
    composition.insert("[UNIMOD:1168]".to_string(), HashMap::from([("H", 0), ("C", 0), ("N", 0), ("O", 2), ("S", 0)])); // Pro->Glu
    composition.insert("[UNIMOD:1169]".to_string(), HashMap::from([("H", 2), ("C", 4), ("N", 0), ("O", 0), ("S", 0)])); // Pro->Phe
    composition.insert("[UNIMOD:1170]".to_string(), HashMap::from([("H", -4), ("C", -3), ("N", 0), ("O", 0), ("S", 0)])); // Pro->Gly
    composition.insert("[UNIMOD:1171]".to_string(), HashMap::from([("H", 5), ("C", 1), ("N", 1), ("O", 0), ("S", 0)])); // Pro->Lys
    composition.insert("[UNIMOD:1172]".to_string(), HashMap::from([("H", 2), ("C", 0), ("N", 0), ("O", 0), ("S", 1)])); // Pro->Met
    composition.insert("[UNIMOD:1173]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", 1), ("O", 1), ("S", 0)])); // Pro->Asn
    composition.insert("[UNIMOD:1174]".to_string(), HashMap::from([("H", 2), ("C", 0), ("N", 0), ("O", 0), ("S", 0)])); // Pro->Val
    composition.insert("[UNIMOD:1175]".to_string(), HashMap::from([("H", 3), ("C", 6), ("N", 1), ("O", 0), ("S", 0)])); // Pro->Trp
    composition.insert("[UNIMOD:1176]".to_string(), HashMap::from([("H", 2), ("C", 4), ("N", 0), ("O", 1), ("S", 0)])); // Pro->Tyr
    composition.insert("[UNIMOD:1177]".to_string(), HashMap::from([("H", -3), ("C", -2), ("N", -1), ("O", -1), ("S", 0)])); // Gln->Ala
    composition.insert("[UNIMOD:1178]".to_string(), HashMap::from([("H", -3), ("C", -2), ("N", -1), ("O", -1), ("S", 1)])); // Gln->Cys
    composition.insert("[UNIMOD:1179]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -1), ("O", 1), ("S", 0)])); // Gln->Asp
    composition.insert("[UNIMOD:1180]".to_string(), HashMap::from([("H", 1), ("C", 4), ("N", -1), ("O", -1), ("S", 0)])); // Gln->Phe
    composition.insert("[UNIMOD:1181]".to_string(), HashMap::from([("H", -5), ("C", -3), ("N", -1), ("O", -1), ("S", 0)])); // Gln->Gly
    composition.insert("[UNIMOD:1182]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", -1), ("O", -1), ("S", 1)])); // Gln->Met
    composition.insert("[UNIMOD:1183]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 0), ("S", 0)])); // Gln->Asn
    composition.insert("[UNIMOD:1184]".to_string(), HashMap::from([("H", -3), ("C", -2), ("N", -1), ("O", 0), ("S", 0)])); // Gln->Ser
    composition.insert("[UNIMOD:1185]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", -1), ("O", 0), ("S", 0)])); // Gln->Thr
    composition.insert("[UNIMOD:1186]".to_string(), HashMap::from([("H", 1), ("C", 0), ("N", -1), ("O", -1), ("S", 0)])); // Gln->Val
    composition.insert("[UNIMOD:1187]".to_string(), HashMap::from([("H", 2), ("C", 6), ("N", 0), ("O", -1), ("S", 0)])); // Gln->Trp
    composition.insert("[UNIMOD:1188]".to_string(), HashMap::from([("H", 1), ("C", 4), ("N", -1), ("O", 0), ("S", 0)])); // Gln->Tyr
    composition.insert("[UNIMOD:1189]".to_string(), HashMap::from([("H", -7), ("C", -3), ("N", -3), ("O", 0), ("S", 0)])); // Arg->Ala
    composition.insert("[UNIMOD:1190]".to_string(), HashMap::from([("H", -7), ("C", -2), ("N", -3), ("O", 2), ("S", 0)])); // Arg->Asp
    composition.insert("[UNIMOD:1191]".to_string(), HashMap::from([("H", -5), ("C", -1), ("N", -3), ("O", 2), ("S", 0)])); // Arg->Glu
    composition.insert("[UNIMOD:1192]".to_string(), HashMap::from([("H", -6), ("C", -2), ("N", -2), ("O", 1), ("S", 0)])); // Arg->Asn
    composition.insert("[UNIMOD:1193]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", -3), ("O", 0), ("S", 0)])); // Arg->Val
    composition.insert("[UNIMOD:1194]".to_string(), HashMap::from([("H", -3), ("C", 3), ("N", -3), ("O", 1), ("S", 0)])); // Arg->Tyr
    composition.insert("[UNIMOD:1195]".to_string(), HashMap::from([("H", -3), ("C", 3), ("N", -3)])); // Arg->Phe
    composition.insert("[UNIMOD:1196]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 0), ("O", 1), ("S", 0)])); // Ser->Asp
    composition.insert("[UNIMOD:1197]".to_string(), HashMap::from([("H", 2), ("C", 2), ("N", 0), ("O", 1), ("S", 0)])); // Ser->Glu
    composition.insert("[UNIMOD:1198]".to_string(), HashMap::from([("H", 2), ("C", 3), ("N", 2), ("O", -1), ("S", 0)])); // Ser->His
    composition.insert("[UNIMOD:1199]".to_string(), HashMap::from([("H", 7), ("C", 3), ("N", 1), ("O", -1), ("S", 0)])); // Ser->Lys
    composition.insert("[UNIMOD:1200]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", -1), ("S", 1)])); // Ser->Met
    composition.insert("[UNIMOD:1201]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 0), ("S", 0)])); // Ser->Gln
    composition.insert("[UNIMOD:1202]".to_string(), HashMap::from([("H", 4), ("C", 2), ("N", 0), ("O", -1), ("S", 0)])); // Ser->Val
    composition.insert("[UNIMOD:1203]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", -1), ("S", 1)])); // Thr->Cys
    composition.insert("[UNIMOD:1204]".to_string(), HashMap::from([("H", -2), ("C", 0), ("N", 0), ("O", 1), ("S", 0)])); // Thr->Asp
    composition.insert("[UNIMOD:1205]".to_string(), HashMap::from([("H", 0), ("C", 1), ("N", 0), ("O", 1), ("S", 0)])); // Thr->Glu
    composition.insert("[UNIMOD:1206]".to_string(), HashMap::from([("H", 2), ("C", 5), ("N", 0), ("O", -1), ("S", 0)])); // Thr->Phe
    composition.insert("[UNIMOD:1207]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", -1), ("S", 0)])); // Thr->Gly
    composition.insert("[UNIMOD:1208]".to_string(), HashMap::from([("H", 0), ("C", 2), ("N", 2), ("O", -1), ("S", 0)])); // Thr->His
    composition.insert("[UNIMOD:1209]".to_string(), HashMap::from([("H", 1), ("C", 1), ("N", 1), ("O", 0), ("S", 0)])); // Thr->Gln
    composition.insert("[UNIMOD:1210]".to_string(), HashMap::from([("H", 2), ("C", 1), ("N", 0), ("O", -1), ("S", 0)])); // Thr->Val
    composition.insert("[UNIMOD:1211]".to_string(), HashMap::from([("H", 3), ("C", 7), ("N", 1), ("O", -1), ("S", 0)])); // Thr->Trp
    composition.insert("[UNIMOD:1212]".to_string(), HashMap::from([("H", 2), ("C", 5), ("N", 0), ("O", 0), ("S", 0)])); // Thr->Tyr
    composition.insert("[UNIMOD:1213]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 0), ("S", 1)])); // Val->Cys
    composition.insert("[UNIMOD:1214]".to_string(), HashMap::from([("H", -2), ("C", 1), ("N", 2), ("O", 0), ("S", 0)])); // Val->His
    composition.insert("[UNIMOD:1215]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 1), ("O", 0), ("S", 0)])); // Val->Lys
    composition.insert("[UNIMOD:1216]".to_string(), HashMap::from([("H", -3), ("C", -1), ("N", 1), ("O", 1), ("S", 0)])); // Val->Asn
    composition.insert("[UNIMOD:1217]".to_string(), HashMap::from([("H", -2), ("C", 0), ("N", 0), ("O", 0), ("S", 0)])); // Val->Pro
    composition.insert("[UNIMOD:1218]".to_string(), HashMap::from([("H", -1), ("C", 0), ("N", 1), ("O", 1), ("S", 0)])); // Val->Gln
    composition.insert("[UNIMOD:1219]".to_string(), HashMap::from([("H", 3), ("C", 1), ("N", 3), ("O", 0), ("S", 0)])); // Val->Arg
    composition.insert("[UNIMOD:1220]".to_string(), HashMap::from([("H", -4), ("C", -2), ("N", 0), ("O", 1), ("S", 0)])); // Val->Ser
    composition.insert("[UNIMOD:1221]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", 0), ("O", 1), ("S", 0)])); // Val->Thr
    composition.insert("[UNIMOD:1222]".to_string(), HashMap::from([("H", 1), ("C", 6), ("N", 1), ("O", 0), ("S", 0)])); // Val->Trp
    composition.insert("[UNIMOD:1223]".to_string(), HashMap::from([("H", 0), ("C", 4), ("N", 0), ("O", 1), ("S", 0)])); // Val->Tyr
    composition.insert("[UNIMOD:1224]".to_string(), HashMap::from([("H", -5), ("C", -8), ("N", -1), ("O", 0), ("S", 0)])); // Trp->Ala
    composition.insert("[UNIMOD:1225]".to_string(), HashMap::from([("H", -5), ("C", -7), ("N", -1), ("O", 2), ("S", 0)])); // Trp->Asp
    composition.insert("[UNIMOD:1226]".to_string(), HashMap::from([("H", -3), ("C", -6), ("N", -1), ("O", 2), ("S", 0)])); // Trp->Glu
    composition.insert("[UNIMOD:1227]".to_string(), HashMap::from([("H", -1), ("C", -2), ("N", -1), ("O", 0), ("S", 0)])); // Trp->Phe
    composition.insert("[UNIMOD:1228]".to_string(), HashMap::from([("H", -3), ("C", -5), ("N", 1), ("O", 0), ("S", 0)])); // Trp->His
    composition.insert("[UNIMOD:1229]".to_string(), HashMap::from([("H", 2), ("C", -5), ("N", 0), ("O", 0), ("S", 0)])); // Trp->Lys
    composition.insert("[UNIMOD:1230]".to_string(), HashMap::from([("H", -1), ("C", -6), ("N", -1), ("O", 0), ("S", 1)])); // Trp->Met
    composition.insert("[UNIMOD:1231]".to_string(), HashMap::from([("H", -4), ("C", -7), ("N", 0), ("O", 1), ("S", 0)])); // Trp->Asn
    composition.insert("[UNIMOD:1232]".to_string(), HashMap::from([("H", -3), ("C", -6), ("N", -1), ("O", 0), ("S", 0)])); // Trp->Pro
    composition.insert("[UNIMOD:1233]".to_string(), HashMap::from([("H", -2), ("C", -6), ("N", 0), ("O", 1), ("S", 0)])); // Trp->Gln
    composition.insert("[UNIMOD:1234]".to_string(), HashMap::from([("H", -3), ("C", -7), ("N", -1), ("O", 1), ("S", 0)])); // Trp->Thr
    composition.insert("[UNIMOD:1235]".to_string(), HashMap::from([("H", -1), ("C", -6), ("N", -1), ("O", 0), ("S", 0)])); // Trp->Val
    composition.insert("[UNIMOD:1236]".to_string(), HashMap::from([("H", -1), ("C", -2), ("N", -1), ("O", 1), ("S", 0)])); // Trp->Tyr
    composition.insert("[UNIMOD:1237]".to_string(), HashMap::from([("H", -4), ("C", -6), ("N", 0), ("O", -1), ("S", 0)])); // Tyr->Ala
    composition.insert("[UNIMOD:1238]".to_string(), HashMap::from([("H", -2), ("C", -4), ("N", 0), ("O", 1), ("S", 0)])); // Tyr->Glu
    composition.insert("[UNIMOD:1239]".to_string(), HashMap::from([("H", -6), ("C", -7), ("N", 0), ("O", -1), ("S", 0)])); // Tyr->Gly
    composition.insert("[UNIMOD:1240]".to_string(), HashMap::from([("H", 3), ("C", -3), ("N", 1), ("O", -1), ("S", 0)])); // Tyr->Lys
    composition.insert("[UNIMOD:1241]".to_string(), HashMap::from([("H", 0), ("C", -4), ("N", 0), ("O", -1), ("S", 1)])); // Tyr->Met
    composition.insert("[UNIMOD:1242]".to_string(), HashMap::from([("H", -2), ("C", -4), ("N", 0), ("O", -1), ("S", 0)])); // Tyr->Pro
    composition.insert("[UNIMOD:1243]".to_string(), HashMap::from([("H", -1), ("C", -4), ("N", 1), ("O", 0), ("S", 0)])); // Tyr->Gln
    composition.insert("[UNIMOD:1244]".to_string(), HashMap::from([("H", 3), ("C", -3), ("N", 3), ("O", -1), ("S", 0)])); // Tyr->Arg
    composition.insert("[UNIMOD:1245]".to_string(), HashMap::from([("H", -2), ("C", -5), ("N", 0), ("O", 0), ("S", 0)])); // Tyr->Thr
    composition.insert("[UNIMOD:1246]".to_string(), HashMap::from([("H", 0), ("C", -4), ("N", 0), ("O", -1), ("S", 0)])); // Tyr->Val
    composition.insert("[UNIMOD:1247]".to_string(), HashMap::from([("H", 1), ("C", 2), ("N", 1), ("O", -1), ("S", 0)])); // Tyr->Trp
    composition.insert("[UNIMOD:1248]".to_string(), HashMap::from([("H", 2), ("C", -3), ("O", -1)])); // Tyr->Xle
    composition.insert("[UNIMOD:1249]".to_string(), HashMap::from([("H", 9), ("C", 7), ("N", 5), ("O", 2)])); // AHA-SS
    composition.insert("[UNIMOD:1250]".to_string(), HashMap::from([("H", 12), ("C", 9), ("N", 6), ("O", 3)])); // AHA-SS_CAM
    composition.insert("[UNIMOD:1251]".to_string(), HashMap::from([("H", 36), ("C", 25), ("N", 6), ("O", 4), ("S", 2)])); // Biotin:Thermo-33033
    composition.insert("[UNIMOD:1252]".to_string(), HashMap::from([("H", 34), ("C", 25), ("N", 6), ("O", 4), ("S", 2)])); // Biotin:Thermo-33033-H
    composition.insert("[UNIMOD:1253]".to_string(), HashMap::from([("H", 6), ("C", 5), ("O", 4)])); // 2-monomethylsuccinyl
    composition.insert("[UNIMOD:1254]".to_string(), HashMap::from([("H", 6), ("C", 7), ("O", 1)])); // Saligenin
    composition.insert("[UNIMOD:1255]".to_string(), HashMap::from([("H", 7), ("C", 7), ("O", 3), ("P", 1)])); // Cresylphosphate
    composition.insert("[UNIMOD:1256]".to_string(), HashMap::from([("H", 13), ("C", 14), ("O", 4), ("P", 1)])); // CresylSaligeninPhosphate
    composition.insert("[UNIMOD:1257]".to_string(), HashMap::from([("H", 8), ("C", 4), ("N", 2), ("O", 1)])); // Ub-Br2
    composition.insert("[UNIMOD:1258]".to_string(), HashMap::from([("H", 12), ("C", 7), ("N", 2), ("O", 3)])); // Ub-VME
    composition.insert("[UNIMOD:1261]".to_string(), HashMap::from([("H", 29), ("C", 31), ("N", 6), ("O", 7)])); // Ub-fluorescein
    composition.insert("[UNIMOD:1262]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 4)])); // 2-dimethylsuccinyl
    composition.insert("[UNIMOD:1263]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 1)])); // Gly
    composition.insert("[UNIMOD:1264]".to_string(), HashMap::from([("H", 13), ("C", 9), ("N", 3), ("O", 5)])); // pupylation
    composition.insert("[UNIMOD:1266]".to_string(), HashMap::from([("C", -4), ("13C", 4)])); // Label:13C(4)
    composition.insert("[UNIMOD:1267]".to_string(), HashMap::from([("C", -4), ("13C", 4), ("O", 1)])); // Label:13C(4)+Oxidation
    composition.insert("[UNIMOD:1270]".to_string(), HashMap::from([("H", 7), ("C", 4), ("N", 1), ("O", 1), ("S", 1)])); // HCysThiolactone
    composition.insert("[UNIMOD:1271]".to_string(), HashMap::from([("H", 7), ("C", 4), ("N", 1), ("O", 2), ("S", 1)])); // HCysteinyl
    composition.insert("[UNIMOD:1276]".to_string(), HashMap::from([("H", 60), ("C", 47), ("N", 23), ("O", 10)])); // UgiJoullie
    composition.insert("[UNIMOD:1277]".to_string(), HashMap::from([("H", 11), ("C", 13), ("N", 3), ("O", 1)])); // Dipyridyl
    composition.insert("[UNIMOD:1278]".to_string(), HashMap::from([("H", 2), ("C", 4), ("O", 1)])); // Furan
    composition.insert("[UNIMOD:1279]".to_string(), HashMap::from([("H", 4), ("C", 8), ("O", 2)])); // Difuran
    composition.insert("[UNIMOD:1281]".to_string(), HashMap::from([("H", 17), ("C", 18), ("N", 1), ("O", 1)])); // BMP-piperidinol
    composition.insert("[UNIMOD:1282]".to_string(), HashMap::from([("H", 10), ("C", 7), ("N", 2), ("O", 2)])); // UgiJoullieProGly
    composition.insert("[UNIMOD:1283]".to_string(), HashMap::from([("H", 20), ("C", 14), ("N", 4), ("O", 4)])); // UgiJoullieProGlyProGly
    composition.insert("[UNIMOD:1286]".to_string(), HashMap::from([("H", 40), ("C", 25), ("N", 2), ("O", 18), ("S", 1)])); // IMEHex(2)NeuAc(1)
    composition.insert("[UNIMOD:1287]".to_string(), HashMap::from([("H", -12), ("C", -6), ("N", -4), ("O", -1)])); // Arg-loss
    composition.insert("[UNIMOD:1288]".to_string(), HashMap::from([("H", 12), ("C", 6), ("N", 4), ("O", 1)])); // Arg
    composition.insert("[UNIMOD:1289]".to_string(), HashMap::from([("H", 6), ("C", 4), ("O", 1)])); // Butyryl
    composition.insert("[UNIMOD:1290]".to_string(), HashMap::from([("H", 6), ("C", 4), ("N", 2), ("O", 2)])); // Dicarbamidomethyl
    composition.insert("[UNIMOD:1291]".to_string(), HashMap::from([("H", -2), ("2H", 6), ("C", 2)])); // Dimethyl:2H(6)
    composition.insert("[UNIMOD:1292]".to_string(), HashMap::from([("H", 14), ("C", 9), ("N", 4), ("O", 4)])); // GGQ
    composition.insert("[UNIMOD:1293]".to_string(), HashMap::from([("H", 21), ("C", 13), ("N", 5), ("O", 6)])); // QTGG
    composition.insert("[UNIMOD:1296]".to_string(), HashMap::from([("C", -3), ("13C", 3)])); // Label:13C(3)
    composition.insert("[UNIMOD:1297]".to_string(), HashMap::from([("C", -3), ("13C", 3), ("N", -1), ("15N", 1)])); // Label:13C(3)15N(1)
    composition.insert("[UNIMOD:1298]".to_string(), HashMap::from([("C", -4), ("13C", 4), ("N", -1), ("15N", 1)])); // Label:13C(4)15N(1)
    composition.insert("[UNIMOD:1299]".to_string(), HashMap::from([("H", -10), ("2H", 10)])); // Label:2H(10)
    composition.insert("[UNIMOD:1300]".to_string(), HashMap::from([("H", -4), ("2H", 4), ("C", -1), ("13C", 1)])); // Label:2H(4)13C(1)
    composition.insert("[UNIMOD:1301]".to_string(), HashMap::from([("H", 12), ("C", 6), ("N", 2), ("O", 1)])); // Lys
    composition.insert("[UNIMOD:1302]".to_string(), HashMap::from([("H", 12), ("C", 1), ("13C", 6), ("15N", 2), ("O", 1)])); // mTRAQ:13C(6)15N(2)
    composition.insert("[UNIMOD:1303]".to_string(), HashMap::from([("H", 17), ("C", 11), ("N", 1), ("O", 8)])); // NeuAc
    composition.insert("[UNIMOD:1304]".to_string(), HashMap::from([("H", 17), ("C", 11), ("N", 1), ("O", 9)])); // NeuGc
    composition.insert("[UNIMOD:1305]".to_string(), HashMap::from([("H", 6), ("C", 3)])); // Propyl
    composition.insert("[UNIMOD:1306]".to_string(), HashMap::from([("2H", 6), ("C", 3)])); // Propyl:2H(6)
    composition.insert("[UNIMOD:1310]".to_string(), HashMap::from([("H", 8), ("C", 9), ("O", 1)])); // Propiophenone
    composition.insert("[UNIMOD:1312]".to_string(), HashMap::from([("H", 6), ("C", 3), ("O", 1)])); // Delta:H(6)C(3)O(1)
    composition.insert("[UNIMOD:1313]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 1)])); // Delta:H(8)C(6)O(1)
    composition.insert("[UNIMOD:1314]".to_string(), HashMap::from([("H", 22), ("C", 13), ("N", 4), ("O", 2), ("S", 1)])); // biotinAcrolein298
    composition.insert("[UNIMOD:1315]".to_string(), HashMap::from([("H", 19), ("C", 18), ("N", 1), ("O", 1)])); // MM-diphenylpentanone
    composition.insert("[UNIMOD:1317]".to_string(), HashMap::from([("H", 18), ("C", 18), ("O", 2)])); // EHD-diphenylpentanone
    composition.insert("[UNIMOD:1320]".to_string(), HashMap::from([("H", 39), ("C", 23), ("N", 5), ("O", 9), ("S", 1)])); // Biotin:Thermo-21901+2H2O
    composition.insert("[UNIMOD:1321]".to_string(), HashMap::from([("H", 15), ("C", 7), ("13C", 1), ("15N", 1), ("18O", 1)])); // DiLeu4plex115
    composition.insert("[UNIMOD:1322]".to_string(), HashMap::from([("H", 13), ("2H", 2), ("C", 8), ("N", 1), ("18O", 1)])); // DiLeu4plex
    composition.insert("[UNIMOD:1323]".to_string(), HashMap::from([("H", 13), ("2H", 2), ("C", 7), ("13C", 1), ("15N", 1), ("O", 1)])); // DiLeu4plex117
    composition.insert("[UNIMOD:1324]".to_string(), HashMap::from([("H", 11), ("2H", 4), ("C", 8), ("N", 1), ("O", 1)])); // DiLeu4plex118
    composition.insert("[UNIMOD:1326]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 1), ("O", 2), ("S", 1)])); // NEMsulfur
    composition.insert("[UNIMOD:1327]".to_string(), HashMap::from([("O", 2), ("S", 1)])); // SulfurDioxide
    composition.insert("[UNIMOD:1328]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 3), ("S", 1)])); // NEMsulfurWater
    composition.insert("[UNIMOD:1330]".to_string(), HashMap::from([("H", 22), ("C", 32), ("N", 2)])); // bisANS-sulfonates
    composition.insert("[UNIMOD:1331]".to_string(), HashMap::from([("H", 2), ("C", 6), ("N", 2), ("O", 4)])); // DNCB_hapten
    composition.insert("[UNIMOD:1340]".to_string(), HashMap::from([("H", 71), ("C", 41), ("N", 5), ("O", 16), ("S", 1)])); // Biotin:Thermo-21911
    composition.insert("[UNIMOD:1341]".to_string(), HashMap::from([("H", 28), ("C", 16), ("N", 4), ("O", 3)])); // iodoTMT
    composition.insert("[UNIMOD:1342]".to_string(), HashMap::from([("H", 28), ("C", 12), ("13C", 4), ("N", 3), ("15N", 1), ("O", 3)])); // iodoTMT6plex
    composition.insert("[UNIMOD:1344]".to_string(), HashMap::from([("H", 11), ("C", 6), ("O", 9), ("P", 1)])); // Phosphogluconoylation
    composition.insert("[UNIMOD:1345]".to_string(), HashMap::from([("H", 4), ("C", 7), ("O", 2)])); // PS_Hapten
    composition.insert("[UNIMOD:1348]".to_string(), HashMap::from([("H", 45), ("C", 37), ("N", 4), ("O", 9), ("S", 2)])); // Cy3-maleimide
    composition.insert("[UNIMOD:1349]".to_string(), HashMap::from([("H", 8), ("C", 8), ("N", 2)])); // benzylguanidine
    composition.insert("[UNIMOD:1350]".to_string(), HashMap::from([("H", 10), ("C", 9), ("N", 2), ("O", 1)])); // CarboxymethylDMAP
    composition.insert("[UNIMOD:1355]".to_string(), HashMap::from([("H", -4), ("O", -1)])); // azole
    composition.insert("[UNIMOD:1356]".to_string(), HashMap::from([("H", 9), ("C", 5), ("O", 7), ("P", 1)])); // phosphoRibosyl
    composition.insert("[UNIMOD:1358]".to_string(), HashMap::from([("H", 4), ("2H", 5), ("C", 6), ("N", 1), ("O", 3)])); // NEM:2H(5)+H2O
    composition.insert("[UNIMOD:1363]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 1)])); // Crotonyl
    composition.insert("[UNIMOD:1364]".to_string(), HashMap::from([("H", 10), ("C", 4), ("N", 1), ("O", 2), ("P", 1)])); // O-Et-N-diMePhospho
    composition.insert("[UNIMOD:1365]".to_string(), HashMap::from([("H", 6), ("C", 2), ("N", 1), ("O", 2), ("P", 1)])); // N-dimethylphosphate
    composition.insert("[UNIMOD:1367]".to_string(), HashMap::from([("H", 20), ("C", 12), ("O", 9)])); // dHex(1)Hex(1)
    composition.insert("[UNIMOD:1368]".to_string(), HashMap::from([("H", -2), ("2H", 6), ("C", 3), ("O", 1)])); // Methyl:2H(3)+Acetyl:2H(3)
    composition.insert("[UNIMOD:1370]".to_string(), HashMap::from([("H", -3), ("2H", 3), ("O", 1)])); // Label:2H(3)+Oxidation
    composition.insert("[UNIMOD:1371]".to_string(), HashMap::from([("H", -3), ("2H", 9), ("C", 3)])); // Trimethyl:2H(9)
    composition.insert("[UNIMOD:1372]".to_string(), HashMap::from([("H", 2), ("13C", 2), ("O", 1)])); // Acetyl:13C(2)
    composition.insert("[UNIMOD:1375]".to_string(), HashMap::from([("H", 30), ("C", 18), ("O", 14)])); // dHex(1)Hex(2)
    composition.insert("[UNIMOD:1376]".to_string(), HashMap::from([("H", 40), ("C", 24), ("O", 19)])); // dHex(1)Hex(3)
    composition.insert("[UNIMOD:1377]".to_string(), HashMap::from([("H", 50), ("C", 30), ("O", 24)])); // dHex(1)Hex(4)
    composition.insert("[UNIMOD:1378]".to_string(), HashMap::from([("H", 60), ("C", 36), ("O", 29)])); // dHex(1)Hex(5)
    composition.insert("[UNIMOD:1379]".to_string(), HashMap::from([("H", 70), ("C", 42), ("O", 34)])); // dHex(1)Hex(6)
    composition.insert("[UNIMOD:1380]".to_string(), HashMap::from([("H", 6), ("C", 3), ("O", 2), ("S", 1)])); // methylsulfonylethyl
    composition.insert("[UNIMOD:1381]".to_string(), HashMap::from([("H", 8), ("C", 4), ("O", 2), ("S", 1)])); // ethylsulfonylethyl
    composition.insert("[UNIMOD:1382]".to_string(), HashMap::from([("H", 8), ("C", 8), ("O", 2), ("S", 1)])); // phenylsulfonylethyl
    composition.insert("[UNIMOD:1383]".to_string(), HashMap::from([("H", 10), ("C", 8), ("N", 1), ("O", 5), ("P", 1)])); // PyridoxalPhosphateH2
    composition.insert("[UNIMOD:1384]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", 3)])); // Homocysteic_acid
    composition.insert("[UNIMOD:1385]".to_string(), HashMap::from([("H", 1), ("N", 1)])); // Hydroxamic_acid
    composition.insert("[UNIMOD:1387]".to_string(), HashMap::from([("H", 5), ("C", 3), ("O", 6), ("P", 1)])); // 3-phosphoglyceryl
    composition.insert("[UNIMOD:1388]".to_string(), HashMap::from([("H", 11), ("C", 5), ("N", 1), ("O", 1)])); // HN2_mustard
    composition.insert("[UNIMOD:1389]".to_string(), HashMap::from([("H", 13), ("C", 6), ("N", 1), ("O", 2)])); // HN3_mustard
    composition.insert("[UNIMOD:1390]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 1), ("O", 3)])); // Oxidation+NEM
    composition.insert("[UNIMOD:1391]".to_string(), HashMap::from([("H", 21), ("C", 27), ("N", 1), ("O", 7)])); // NHS-fluorescein
    composition.insert("[UNIMOD:1392]".to_string(), HashMap::from([("H", 20), ("C", 7), ("13C", 4), ("N", 1), ("15N", 1), ("O", 2)])); // DiART6plex
    composition.insert("[UNIMOD:1393]".to_string(), HashMap::from([("H", 20), ("C", 8), ("13C", 3), ("15N", 2), ("O", 2)])); // DiART6plex115
    composition.insert("[UNIMOD:1394]".to_string(), HashMap::from([("H", 18), ("2H", 2), ("C", 9), ("13C", 2), ("N", 1), ("15N", 1), ("O", 2)])); // DiART6plex116/119
    composition.insert("[UNIMOD:1395]".to_string(), HashMap::from([("H", 18), ("2H", 2), ("C", 10), ("13C", 1), ("15N", 2), ("O", 2)])); // DiART6plex117
    composition.insert("[UNIMOD:1396]".to_string(), HashMap::from([("H", 18), ("2H", 2), ("C", 8), ("13C", 3), ("N", 2), ("O", 2)])); // DiART6plex118
    composition.insert("[UNIMOD:1397]".to_string(), HashMap::from([("H", 7), ("C", 8), ("N", 1), ("O", 1)])); // Iodoacetanilide
    composition.insert("[UNIMOD:1398]".to_string(), HashMap::from([("H", 7), ("C", 2), ("13C", 6), ("N", 1), ("O", 1)])); // Iodoacetanilide:13C(6)
    composition.insert("[UNIMOD:1399]".to_string(), HashMap::from([("H", 20), ("C", 13), ("N", 2), ("O", 6), ("S", 2)])); // Dap-DSP
    composition.insert("[UNIMOD:1400]".to_string(), HashMap::from([("H", 17), ("C", 11), ("N", 1), ("O", 7)])); // MurNAc
    composition.insert("[UNIMOD:1402]".to_string(), HashMap::from([("H", -7), ("2H", 7), ("N", -4), ("15N", 4)])); // Label:2H(7)15N(4)
    composition.insert("[UNIMOD:1403]".to_string(), HashMap::from([("H", -6), ("2H", 6), ("N", -1), ("15N", 1)])); // Label:2H(6)15N(1)
    composition.insert("[UNIMOD:1405]".to_string(), HashMap::from([("H", 107), ("C", 72), ("N", 17), ("O", 31)])); // EEEDVIEVYQEQTGG
    composition.insert("[UNIMOD:1406]".to_string(), HashMap::from([("H", 102), ("C", 69), ("N", 18), ("O", 30)])); // EDEDTIDVFQQQTGG
    composition.insert("[UNIMOD:1408]".to_string(), HashMap::from([("H", 136), ("C", 84), ("N", 6), ("O", 61)])); // Hex(5)HexNAc(4)NeuAc(2)
    composition.insert("[UNIMOD:1409]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 53)])); // Hex(5)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1410]".to_string(), HashMap::from([("H", 129), ("C", 79), ("N", 5), ("O", 57)])); // dHex(1)Hex(5)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1411]".to_string(), HashMap::from([("H", 146), ("C", 90), ("N", 6), ("O", 65)])); // dHex(1)Hex(5)HexNAc(4)NeuAc(2)
    composition.insert("[UNIMOD:1412]".to_string(), HashMap::from([("H", 13), ("C", 8), ("N", 1), ("O", 8), ("S", 1)])); // s-GlcNAc
    composition.insert("[UNIMOD:1413]".to_string(), HashMap::from([("H", 21), ("C", 12), ("O", 13), ("P", 1)])); // PhosphoHex(2)
    composition.insert("[UNIMOD:1414]".to_string(), HashMap::from([("H", -3), ("2H", 9), ("13C", 3)])); // Trimethyl:13C(3)2H(9)
    composition.insert("[UNIMOD:1419]".to_string(), HashMap::from([("H", -3), ("15N", -1)])); // 15N-oxobutanoic
    composition.insert("[UNIMOD:1420]".to_string(), HashMap::from([("H", 23), ("C", 10), ("N", 3)])); // spermine
    composition.insert("[UNIMOD:1421]".to_string(), HashMap::from([("H", 16), ("C", 7), ("N", 2)])); // spermidine
    composition.insert("[UNIMOD:1423]".to_string(), HashMap::from([("H", 35), ("C", 21), ("N", 3), ("O", 7), ("S", 1)])); // Biotin:Thermo-21330
    composition.insert("[UNIMOD:1425]".to_string(), HashMap::from([("H", 8), ("C", 5), ("O", 4)])); // Pentose
    composition.insert("[UNIMOD:1426]".to_string(), HashMap::from([("H", 18), ("C", 11), ("O", 9)])); // Hex(1)Pent(1)
    composition.insert("[UNIMOD:1427]".to_string(), HashMap::from([("H", 18), ("C", 12), ("O", 11)])); // Hex(1)HexA(1)
    composition.insert("[UNIMOD:1428]".to_string(), HashMap::from([("H", 26), ("C", 16), ("O", 13)])); // Hex(1)Pent(2)
    composition.insert("[UNIMOD:1429]".to_string(), HashMap::from([("H", 24), ("C", 14), ("N", 1), ("O", 13), ("P", 1)])); // Hex(1)HexNAc(1)Phos(1)
    composition.insert("[UNIMOD:1430]".to_string(), HashMap::from([("H", 23), ("C", 14), ("N", 1), ("O", 13), ("S", 1)])); // Hex(1)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1431]".to_string(), HashMap::from([("H", 27), ("C", 17), ("N", 1), ("O", 13)])); // Hex(1)NeuAc(1)
    composition.insert("[UNIMOD:1432]".to_string(), HashMap::from([("H", 27), ("C", 17), ("N", 1), ("O", 14)])); // Hex(1)NeuGc(1)
    composition.insert("[UNIMOD:1433]".to_string(), HashMap::from([("H", 39), ("C", 24), ("N", 3), ("O", 15)])); // HexNAc(3)
    composition.insert("[UNIMOD:1434]".to_string(), HashMap::from([("H", 30), ("C", 19), ("N", 2), ("O", 13)])); // HexNAc(1)NeuAc(1)
    composition.insert("[UNIMOD:1435]".to_string(), HashMap::from([("H", 30), ("C", 19), ("N", 2), ("O", 14)])); // HexNAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1436]".to_string(), HashMap::from([("H", 35), ("C", 21), ("N", 1), ("O", 14)])); // Hex(1)HexNAc(1)dHex(1)Me(1)
    composition.insert("[UNIMOD:1437]".to_string(), HashMap::from([("H", 37), ("C", 22), ("N", 1), ("O", 14)])); // Hex(1)HexNAc(1)dHex(1)Me(2)
    composition.insert("[UNIMOD:1438]".to_string(), HashMap::from([("H", 33), ("C", 20), ("N", 1), ("O", 15)])); // Hex(2)HexNAc(1)
    composition.insert("[UNIMOD:1439]".to_string(), HashMap::from([("H", 31), ("C", 20), ("N", 1), ("O", 16)])); // Hex(1)HexA(1)HexNAc(1)
    composition.insert("[UNIMOD:1440]".to_string(), HashMap::from([("H", 35), ("C", 21), ("N", 1), ("O", 15)])); // Hex(2)HexNAc(1)Me(1)
    composition.insert("[UNIMOD:1441]".to_string(), HashMap::from([("H", 34), ("C", 21), ("O", 17)])); // Hex(1)Pent(3)
    composition.insert("[UNIMOD:1442]".to_string(), HashMap::from([("H", 35), ("C", 22), ("N", 1), ("O", 17)])); // Hex(1)NeuAc(1)Pent(1)
    composition.insert("[UNIMOD:1443]".to_string(), HashMap::from([("H", 33), ("C", 20), ("N", 1), ("O", 18), ("S", 1)])); // Hex(2)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1444]".to_string(), HashMap::from([("H", 37), ("C", 23), ("N", 1), ("O", 18)])); // Hex(2)NeuAc(1)
    composition.insert("[UNIMOD:1445]".to_string(), HashMap::from([("H", 40), ("C", 24), ("O", 18)])); // dHex(2)Hex(2)
    composition.insert("[UNIMOD:1446]".to_string(), HashMap::from([("H", 38), ("C", 24), ("O", 20)])); // dHex(1)Hex(2)HexA(1)
    composition.insert("[UNIMOD:1447]".to_string(), HashMap::from([("H", 36), ("C", 22), ("N", 2), ("O", 18), ("S", 1)])); // Hex(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1448]".to_string(), HashMap::from([("H", 40), ("C", 24), ("O", 20)])); // Hex(4)
    composition.insert("[UNIMOD:1449]".to_string(), HashMap::from([("H", 64), ("C", 39), ("N", 2), ("O", 28)])); // dHex(1)Hex(2)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1450]".to_string(), HashMap::from([("H", 63), ("C", 39), ("N", 3), ("O", 28)])); // Hex(2)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1451]".to_string(), HashMap::from([("H", 64), ("C", 39), ("N", 2), ("O", 29)])); // Hex(3)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1452]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 30)])); // Hex(4)HexNAc(2)
    composition.insert("[UNIMOD:1453]".to_string(), HashMap::from([("H", 71), ("C", 43), ("N", 1), ("O", 33)])); // dHex(1)Hex(4)HexNAc(1)Pent(1)
    composition.insert("[UNIMOD:1454]".to_string(), HashMap::from([("H", 74), ("C", 45), ("N", 2), ("O", 33)])); // dHex(1)Hex(3)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1455]".to_string(), HashMap::from([("H", 73), ("C", 45), ("N", 3), ("O", 33)])); // Hex(3)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1456]".to_string(), HashMap::from([("H", 74), ("C", 45), ("N", 2), ("O", 34)])); // Hex(4)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1457]".to_string(), HashMap::from([("H", 77), ("C", 47), ("N", 3), ("O", 34)])); // Hex(3)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1458]".to_string(), HashMap::from([("H", 77), ("C", 46), ("N", 2), ("O", 38), ("P", 1)])); // Hex(5)HexNAc(2)Phos(1)
    composition.insert("[UNIMOD:1459]".to_string(), HashMap::from([("H", 84), ("C", 51), ("N", 2), ("O", 38)])); // dHex(1)Hex(4)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1460]".to_string(), HashMap::from([("H", 83), ("C", 50), ("N", 1), ("O", 40)])); // Hex(7)HexNAc(1)
    composition.insert("[UNIMOD:1461]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 38)])); // Hex(4)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1462]".to_string(), HashMap::from([("H", 86), ("C", 52), ("N", 2), ("O", 39)])); // dHex(1)Hex(5)HexNAc(2)
    composition.insert("[UNIMOD:1463]".to_string(), HashMap::from([("H", 87), ("C", 53), ("N", 3), ("O", 38)])); // dHex(1)Hex(3)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1464]".to_string(), HashMap::from([("H", 82), ("C", 50), ("N", 4), ("O", 38), ("S", 1)])); // Hex(3)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1465]".to_string(), HashMap::from([("H", 86), ("C", 52), ("N", 2), ("O", 40)])); // Hex(6)HexNAc(2)
    composition.insert("[UNIMOD:1466]".to_string(), HashMap::from([("H", 87), ("C", 53), ("N", 3), ("O", 39)])); // Hex(4)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1467]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 39)])); // dHex(1)Hex(4)HexNAc(3)
    composition.insert("[UNIMOD:1468]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 40)])); // Hex(5)HexNAc(3)
    composition.insert("[UNIMOD:1469]".to_string(), HashMap::from([("H", 90), ("C", 55), ("N", 4), ("O", 39)])); // Hex(3)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1470]".to_string(), HashMap::from([("H", 87), ("C", 52), ("N", 2), ("O", 43), ("P", 1)])); // Hex(6)HexNAc(2)Phos(1)
    composition.insert("[UNIMOD:1471]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 42), ("S", 1)])); // dHex(1)Hex(4)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1472]".to_string(), HashMap::from([("H", 94), ("C", 57), ("N", 2), ("O", 43)])); // dHex(1)Hex(5)HexNAc(2)Pent(1)
    composition.insert("[UNIMOD:1473]".to_string(), HashMap::from([("H", 93), ("C", 56), ("N", 1), ("O", 45)])); // Hex(8)HexNAc(1)
    composition.insert("[UNIMOD:1474]".to_string(), HashMap::from([("H", 95), ("C", 58), ("N", 3), ("O", 42)])); // dHex(1)Hex(3)HexNAc(3)Pent(2)
    composition.insert("[UNIMOD:1475]".to_string(), HashMap::from([("H", 97), ("C", 59), ("N", 3), ("O", 42)])); // dHex(2)Hex(3)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1476]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 42), ("S", 1)])); // dHex(1)Hex(3)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1477]".to_string(), HashMap::from([("H", 96), ("C", 58), ("N", 2), ("O", 44)])); // dHex(1)Hex(6)HexNAc(2)
    composition.insert("[UNIMOD:1478]".to_string(), HashMap::from([("H", 97), ("C", 59), ("N", 3), ("O", 43)])); // dHex(1)Hex(4)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1479]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 43), ("S", 1)])); // Hex(4)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1480]".to_string(), HashMap::from([("H", 96), ("C", 58), ("N", 2), ("O", 45)])); // Hex(7)HexNAc(2)
    composition.insert("[UNIMOD:1481]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 43)])); // dHex(2)Hex(4)HexNAc(3)
    composition.insert("[UNIMOD:1482]".to_string(), HashMap::from([("H", 97), ("C", 59), ("N", 3), ("O", 44)])); // Hex(5)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1483]".to_string(), HashMap::from([("H", 96), ("C", 59), ("N", 4), ("O", 44)])); // Hex(4)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1484]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 44)])); // dHex(1)Hex(5)HexNAc(3)
    composition.insert("[UNIMOD:1485]".to_string(), HashMap::from([("H", 100), ("C", 61), ("N", 4), ("O", 43)])); // dHex(1)Hex(3)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1486]".to_string(), HashMap::from([("H", 95), ("C", 58), ("N", 5), ("O", 43), ("S", 1)])); // Hex(3)HexNAc(5)Sulf(1)
    composition.insert("[UNIMOD:1487]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 45)])); // Hex(6)HexNAc(3)
    composition.insert("[UNIMOD:1488]".to_string(), HashMap::from([("H", 99), ("C", 61), ("N", 5), ("O", 43)])); // Hex(3)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1489]".to_string(), HashMap::from([("H", 100), ("C", 61), ("N", 4), ("O", 44)])); // Hex(4)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1490]".to_string(), HashMap::from([("H", 97), ("C", 58), ("N", 2), ("O", 48), ("P", 1)])); // Hex(7)HexNAc(2)Phos(1)
    composition.insert("[UNIMOD:1491]".to_string(), HashMap::from([("H", 104), ("C", 63), ("N", 4), ("O", 44)])); // Hex(4)HexNAc(4)Me(2)Pent(1)
    composition.insert("[UNIMOD:1492]".to_string(), HashMap::from([("H", 103), ("C", 63), ("N", 3), ("O", 46)])); // dHex(1)Hex(3)HexNAc(3)Pent(3)
    composition.insert("[UNIMOD:1493]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 47), ("S", 1)])); // dHex(1)Hex(5)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1494]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 3), ("O", 46)])); // dHex(2)Hex(3)HexNAc(3)Pent(2)
    composition.insert("[UNIMOD:1495]".to_string(), HashMap::from([("H", 100), ("C", 60), ("N", 3), ("O", 48), ("P", 1)])); // Hex(6)HexNAc(3)Phos(1)
    composition.insert("[UNIMOD:1496]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 5), ("O", 45)])); // Hex(4)HexNAc(5)
    composition.insert("[UNIMOD:1497]".to_string(), HashMap::from([("H", 107), ("C", 65), ("N", 3), ("O", 46)])); // dHex(3)Hex(3)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1498]".to_string(), HashMap::from([("H", 107), ("C", 65), ("N", 3), ("O", 47)])); // dHex(2)Hex(4)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1499]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 47), ("S", 1)])); // dHex(1)Hex(4)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1500]".to_string(), HashMap::from([("H", 106), ("C", 64), ("N", 2), ("O", 49)])); // dHex(1)Hex(7)HexNAc(2)
    composition.insert("[UNIMOD:1501]".to_string(), HashMap::from([("H", 106), ("C", 65), ("N", 4), ("O", 47)])); // dHex(1)Hex(4)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1502]".to_string(), HashMap::from([("H", 98), ("C", 58), ("N", 2), ("O", 51), ("P", 2)])); // Hex(7)HexNAc(2)Phos(2)
    composition.insert("[UNIMOD:1503]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 48), ("S", 1)])); // Hex(5)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1504]".to_string(), HashMap::from([("H", 106), ("C", 64), ("N", 2), ("O", 50)])); // Hex(8)HexNAc(2)
    composition.insert("[UNIMOD:1505]".to_string(), HashMap::from([("H", 108), ("C", 66), ("N", 4), ("O", 47)])); // dHex(1)Hex(3)HexNAc(4)Pent(2)
    composition.insert("[UNIMOD:1506]".to_string(), HashMap::from([("H", 106), ("C", 65), ("N", 4), ("O", 48)])); // dHex(1)Hex(4)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1507]".to_string(), HashMap::from([("H", 110), ("C", 67), ("N", 4), ("O", 47)])); // dHex(2)Hex(3)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1508]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 5), ("O", 47), ("S", 1)])); // dHex(1)Hex(3)HexNAc(5)Sulf(1)
    composition.insert("[UNIMOD:1509]".to_string(), HashMap::from([("H", 109), ("C", 66), ("N", 3), ("O", 49)])); // dHex(1)Hex(6)HexNAc(3)
    composition.insert("[UNIMOD:1510]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 47)])); // dHex(1)Hex(3)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1511]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 47)])); // dHex(3)Hex(3)HexNAc(4)
    composition.insert("[UNIMOD:1512]".to_string(), HashMap::from([("H", 110), ("C", 67), ("N", 4), ("O", 48)])); // dHex(1)Hex(4)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1513]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 5), ("O", 48), ("S", 1)])); // Hex(4)HexNAc(5)Sulf(1)
    composition.insert("[UNIMOD:1514]".to_string(), HashMap::from([("H", 109), ("C", 66), ("N", 3), ("O", 50)])); // Hex(7)HexNAc(3)
    composition.insert("[UNIMOD:1515]".to_string(), HashMap::from([("H", 106), ("C", 65), ("N", 4), ("O", 50), ("S", 1)])); // dHex(1)Hex(4)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1516]".to_string(), HashMap::from([("H", 114), ("C", 69), ("N", 4), ("O", 49)])); // Hex(5)HexNAc(4)Me(2)Pent(1)
    composition.insert("[UNIMOD:1517]".to_string(), HashMap::from([("H", 108), ("C", 66), ("N", 6), ("O", 48), ("S", 1)])); // Hex(3)HexNAc(6)Sulf(1)
    composition.insert("[UNIMOD:1518]".to_string(), HashMap::from([("H", 109), ("C", 66), ("N", 3), ("O", 52), ("S", 1)])); // dHex(1)Hex(6)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1519]".to_string(), HashMap::from([("H", 115), ("C", 70), ("N", 5), ("O", 49)])); // dHex(1)Hex(4)HexNAc(5)
    composition.insert("[UNIMOD:1520]".to_string(), HashMap::from([("H", 107), ("C", 66), ("N", 3), ("O", 53), ("S", 1)])); // dHex(1)Hex(5)HexA(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1521]".to_string(), HashMap::from([("H", 110), ("C", 66), ("N", 3), ("O", 53), ("P", 1)])); // Hex(7)HexNAc(3)Phos(1)
    composition.insert("[UNIMOD:1522]".to_string(), HashMap::from([("H", 118), ("C", 71), ("N", 4), ("O", 50)])); // Hex(6)HexNAc(4)Me(3)
    composition.insert("[UNIMOD:1523]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 51), ("S", 1)])); // dHex(2)Hex(4)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1524]".to_string(), HashMap::from([("H", 113), ("C", 70), ("N", 5), ("O", 51)])); // Hex(4)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1525]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 51)])); // dHex(1)Hex(3)HexNAc(4)Pent(3)
    composition.insert("[UNIMOD:1526]".to_string(), HashMap::from([("H", 117), ("C", 71), ("N", 3), ("O", 52)])); // dHex(2)Hex(5)HexNAc(3)Pent(1)
    composition.insert("[UNIMOD:1527]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 52), ("S", 1)])); // dHex(1)Hex(5)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1528]".to_string(), HashMap::from([("H", 118), ("C", 72), ("N", 4), ("O", 51)])); // dHex(2)Hex(3)HexNAc(4)Pent(2)
    composition.insert("[UNIMOD:1529]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 52)])); // dHex(1)Hex(5)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1530]".to_string(), HashMap::from([("H", 108), ("C", 66), ("N", 6), ("O", 51), ("S", 2)])); // Hex(3)HexNAc(6)Sulf(2)
    composition.insert("[UNIMOD:1531]".to_string(), HashMap::from([("H", 116), ("C", 70), ("N", 2), ("O", 55)])); // Hex(9)HexNAc(2)
    composition.insert("[UNIMOD:1532]".to_string(), HashMap::from([("H", 118), ("C", 72), ("N", 6), ("O", 50)])); // Hex(4)HexNAc(6)
    composition.insert("[UNIMOD:1533]".to_string(), HashMap::from([("H", 120), ("C", 73), ("N", 4), ("O", 51)])); // dHex(3)Hex(3)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1534]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 53)])); // dHex(1)Hex(5)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1535]".to_string(), HashMap::from([("H", 120), ("C", 73), ("N", 4), ("O", 52)])); // dHex(2)Hex(4)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1536]".to_string(), HashMap::from([("H", 115), ("C", 70), ("N", 5), ("O", 52), ("S", 1)])); // dHex(1)Hex(4)HexNAc(5)Sulf(1)
    composition.insert("[UNIMOD:1537]".to_string(), HashMap::from([("H", 119), ("C", 72), ("N", 3), ("O", 54)])); // dHex(1)Hex(7)HexNAc(3)
    composition.insert("[UNIMOD:1538]".to_string(), HashMap::from([("H", 120), ("C", 73), ("N", 4), ("O", 53)])); // dHex(1)Hex(5)HexNAc(4)Pent(1)
    composition.insert("[UNIMOD:1539]".to_string(), HashMap::from([("H", 107), ("C", 66), ("N", 3), ("O", 56), ("S", 2)])); // dHex(1)Hex(5)HexA(1)HexNAc(3)Sulf(2)
    composition.insert("[UNIMOD:1540]".to_string(), HashMap::from([("H", 121), ("C", 74), ("N", 7), ("O", 50)])); // Hex(3)HexNAc(7)
    composition.insert("[UNIMOD:1541]".to_string(), HashMap::from([("H", 122), ("C", 74), ("N", 4), ("O", 53)])); // dHex(2)Hex(5)HexNAc(4)
    composition.insert("[UNIMOD:1542]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 54), ("S", 1)])); // dHex(2)Hex(4)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1543]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 55), ("S", 2)])); // dHex(1)Hex(5)HexNAc(4)Sulf(2)
    composition.insert("[UNIMOD:1544]".to_string(), HashMap::from([("H", 124), ("C", 75), ("N", 4), ("O", 53)])); // dHex(1)Hex(5)HexNAc(4)Me(2)Pent(1)
    composition.insert("[UNIMOD:1545]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 54)])); // Hex(5)HexNAc(4)NeuGc(1)
    composition.insert("[UNIMOD:1546]".to_string(), HashMap::from([("H", 118), ("C", 72), ("N", 6), ("O", 52), ("S", 1)])); // dHex(1)Hex(3)HexNAc(6)Sulf(1)
    composition.insert("[UNIMOD:1547]".to_string(), HashMap::from([("H", 122), ("C", 74), ("N", 4), ("O", 54)])); // dHex(1)Hex(6)HexNAc(4)
    composition.insert("[UNIMOD:1548]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 55), ("S", 1)])); // dHex(1)Hex(5)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1549]".to_string(), HashMap::from([("H", 122), ("C", 74), ("N", 4), ("O", 55)])); // Hex(7)HexNAc(4)
    composition.insert("[UNIMOD:1550]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 56), ("S", 1)])); // dHex(1)Hex(5)HexNAc(3)NeuGc(1)Sulf(1)
    composition.insert("[UNIMOD:1551]".to_string(), HashMap::from([("H", 122), ("C", 75), ("N", 6), ("O", 53)])); // Hex(4)HexNAc(5)NeuAc(1)
    composition.insert("[UNIMOD:1552]".to_string(), HashMap::from([("H", 126), ("C", 76), ("N", 4), ("O", 54)])); // Hex(6)HexNAc(4)Me(3)Pent(1)
    composition.insert("[UNIMOD:1553]".to_string(), HashMap::from([("H", 119), ("C", 72), ("N", 3), ("O", 57), ("S", 1)])); // dHex(1)Hex(7)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1554]".to_string(), HashMap::from([("H", 120), ("C", 72), ("N", 3), ("O", 57), ("P", 1)])); // dHex(1)Hex(7)HexNAc(3)Phos(1)
    composition.insert("[UNIMOD:1555]".to_string(), HashMap::from([("H", 125), ("C", 76), ("N", 5), ("O", 54)])); // dHex(1)Hex(5)HexNAc(5)
    composition.insert("[UNIMOD:1556]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 55), ("S", 1)])); // dHex(1)Hex(4)HexNAc(4)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1557]".to_string(), HashMap::from([("H", 122), ("C", 74), ("N", 4), ("O", 55), ("S", 1)])); // dHex(3)Hex(4)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1558]".to_string(), HashMap::from([("H", 121), ("C", 74), ("N", 7), ("O", 53), ("S", 1)])); // Hex(3)HexNAc(7)Sulf(1)
    composition.insert("[UNIMOD:1559]".to_string(), HashMap::from([("H", 125), ("C", 76), ("N", 5), ("O", 55)])); // Hex(6)HexNAc(5)
    composition.insert("[UNIMOD:1560]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 56), ("S", 1)])); // Hex(5)HexNAc(4)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1561]".to_string(), HashMap::from([("H", 125), ("C", 77), ("N", 7), ("O", 53)])); // Hex(3)HexNAc(6)NeuAc(1)
    composition.insert("[UNIMOD:1562]".to_string(), HashMap::from([("H", 128), ("C", 78), ("N", 6), ("O", 53)])); // dHex(2)Hex(3)HexNAc(6)
    composition.insert("[UNIMOD:1563]".to_string(), HashMap::from([("H", 40), ("C", 25), ("N", 2), ("O", 19)])); // Hex(1)HexNAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1564]".to_string(), HashMap::from([("H", 43), ("C", 26), ("N", 1), ("O", 19)])); // dHex(1)Hex(2)HexNAc(1)
    composition.insert("[UNIMOD:1565]".to_string(), HashMap::from([("H", 39), ("C", 24), ("N", 3), ("O", 18), ("S", 1)])); // HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1566]".to_string(), HashMap::from([("H", 43), ("C", 26), ("N", 1), ("O", 20)])); // Hex(3)HexNAc(1)
    composition.insert("[UNIMOD:1567]".to_string(), HashMap::from([("H", 37), ("C", 23), ("N", 1), ("O", 21), ("S", 1)])); // Hex(1)HexNAc(1)Kdn(1)Sulf(1)
    composition.insert("[UNIMOD:1568]".to_string(), HashMap::from([("H", 43), ("C", 27), ("N", 3), ("O", 18)])); // HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1570]".to_string(), HashMap::from([("H", 41), ("C", 26), ("N", 1), ("O", 21)])); // HexNAc(1)Kdn(2)
    composition.insert("[UNIMOD:1571]".to_string(), HashMap::from([("H", 45), ("C", 27), ("N", 1), ("O", 20)])); // Hex(3)HexNAc(1)Me(1)
    composition.insert("[UNIMOD:1572]".to_string(), HashMap::from([("H", 36), ("C", 23), ("O", 23), ("S", 1)])); // Hex(2)HexA(1)Pent(1)Sulf(1)
    composition.insert("[UNIMOD:1573]".to_string(), HashMap::from([("H", 43), ("C", 27), ("N", 3), ("O", 19)])); // HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1575]".to_string(), HashMap::from([("H", 41), ("C", 24), ("O", 23), ("P", 1)])); // Hex(4)Phos(1)
    composition.insert("[UNIMOD:1577]".to_string(), HashMap::from([("H", 40), ("C", 25), ("N", 2), ("O", 21), ("S", 1)])); // Hex(1)HexNAc(1)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1578]".to_string(), HashMap::from([("H", 44), ("C", 28), ("N", 2), ("O", 21)])); // Hex(1)HexA(1)HexNAc(2)
    composition.insert("[UNIMOD:1579]".to_string(), HashMap::from([("H", 43), ("C", 26), ("N", 1), ("O", 22), ("S", 1)])); // dHex(1)Hex(2)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1580]".to_string(), HashMap::from([("H", 49), ("C", 30), ("N", 3), ("O", 19)])); // dHex(1)HexNAc(3)
    composition.insert("[UNIMOD:1581]".to_string(), HashMap::from([("H", 47), ("C", 29), ("N", 1), ("O", 22)])); // dHex(1)Hex(1)HexNAc(1)Kdn(1)
    composition.insert("[UNIMOD:1582]".to_string(), HashMap::from([("H", 49), ("C", 30), ("N", 3), ("O", 20)])); // Hex(1)HexNAc(3)
    composition.insert("[UNIMOD:1583]".to_string(), HashMap::from([("H", 43), ("C", 27), ("N", 3), ("O", 21), ("S", 1)])); // HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1584]".to_string(), HashMap::from([("H", 50), ("C", 30), ("O", 23)])); // dHex(2)Hex(3)
    composition.insert("[UNIMOD:1585]".to_string(), HashMap::from([("H", 41), ("C", 26), ("N", 1), ("O", 24), ("S", 1)])); // Hex(2)HexA(1)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1586]".to_string(), HashMap::from([("H", 48), ("C", 30), ("O", 24)])); // dHex(2)Hex(2)HexA(1)
    composition.insert("[UNIMOD:1587]".to_string(), HashMap::from([("H", 46), ("C", 28), ("N", 2), ("O", 22), ("S", 1)])); // dHex(1)Hex(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1588]".to_string(), HashMap::from([("H", 50), ("C", 31), ("N", 2), ("O", 22)])); // dHex(1)Hex(1)HexNAc(1)NeuAc(1)
    composition.insert("[UNIMOD:1589]".to_string(), HashMap::from([("H", 46), ("C", 28), ("N", 2), ("O", 23), ("S", 1)])); // Hex(2)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1590]".to_string(), HashMap::from([("H", 50), ("C", 30), ("O", 25)])); // Hex(5)
    composition.insert("[UNIMOD:1591]".to_string(), HashMap::from([("H", 52), ("C", 32), ("N", 4), ("O", 20)])); // HexNAc(4)
    composition.insert("[UNIMOD:1592]".to_string(), HashMap::from([("H", 47), ("C", 30), ("N", 3), ("O", 23)])); // HexNAc(1)NeuGc(2)
    composition.insert("[UNIMOD:1593]".to_string(), HashMap::from([("H", 50), ("C", 31), ("N", 2), ("O", 23)])); // dHex(1)Hex(1)HexNAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1594]".to_string(), HashMap::from([("H", 53), ("C", 32), ("N", 1), ("O", 23)])); // dHex(2)Hex(2)HexNAc(1)
    composition.insert("[UNIMOD:1595]".to_string(), HashMap::from([("H", 50), ("C", 31), ("N", 2), ("O", 24)])); // Hex(2)HexNAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1596]".to_string(), HashMap::from([("H", 53), ("C", 32), ("N", 1), ("O", 24)])); // dHex(1)Hex(3)HexNAc(1)
    composition.insert("[UNIMOD:1597]".to_string(), HashMap::from([("H", 51), ("C", 32), ("N", 1), ("O", 25)])); // dHex(1)Hex(2)HexA(1)HexNAc(1)
    composition.insert("[UNIMOD:1598]".to_string(), HashMap::from([("H", 49), ("C", 30), ("N", 3), ("O", 23), ("S", 1)])); // Hex(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1599]".to_string(), HashMap::from([("H", 53), ("C", 32), ("N", 1), ("O", 25)])); // Hex(4)HexNAc(1)
    composition.insert("[UNIMOD:1600]".to_string(), HashMap::from([("H", 53), ("C", 33), ("N", 3), ("O", 23)])); // Hex(1)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1602]".to_string(), HashMap::from([("H", 53), ("C", 33), ("N", 3), ("O", 24)])); // Hex(1)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1604]".to_string(), HashMap::from([("H", 51), ("C", 30), ("O", 28), ("P", 1)])); // Hex(5)Phos(1)
    composition.insert("[UNIMOD:1606]".to_string(), HashMap::from([("H", 57), ("C", 35), ("N", 1), ("O", 26)])); // dHex(2)Hex(1)HexNAc(1)Kdn(1)
    composition.insert("[UNIMOD:1607]".to_string(), HashMap::from([("H", 53), ("C", 32), ("N", 1), ("O", 27), ("S", 1)])); // dHex(1)Hex(3)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1608]".to_string(), HashMap::from([("H", 59), ("C", 36), ("N", 3), ("O", 24)])); // dHex(1)Hex(1)HexNAc(3)
    composition.insert("[UNIMOD:1609]".to_string(), HashMap::from([("H", 51), ("C", 32), ("N", 1), ("O", 28), ("S", 1)])); // dHex(1)Hex(2)HexA(1)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1610]".to_string(), HashMap::from([("H", 59), ("C", 36), ("N", 3), ("O", 25)])); // Hex(2)HexNAc(3)
    composition.insert("[UNIMOD:1611]".to_string(), HashMap::from([("H", 53), ("C", 33), ("N", 3), ("O", 26), ("S", 1)])); // Hex(1)HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1612]".to_string(), HashMap::from([("H", 60), ("C", 36), ("O", 28)])); // dHex(2)Hex(4)
    composition.insert("[UNIMOD:1614]".to_string(), HashMap::from([("H", 60), ("C", 37), ("N", 2), ("O", 26)])); // dHex(2)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1615]".to_string(), HashMap::from([("H", 56), ("C", 34), ("N", 2), ("O", 27), ("S", 1)])); // dHex(1)Hex(2)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1616]".to_string(), HashMap::from([("H", 62), ("C", 38), ("N", 4), ("O", 24)])); // dHex(1)HexNAc(4)
    composition.insert("[UNIMOD:1617]".to_string(), HashMap::from([("H", 57), ("C", 36), ("N", 3), ("O", 27)])); // Hex(1)HexNAc(1)NeuAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1618]".to_string(), HashMap::from([("H", 60), ("C", 37), ("N", 2), ("O", 27)])); // dHex(1)Hex(1)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1619]".to_string(), HashMap::from([("H", 57), ("C", 36), ("N", 3), ("O", 28)])); // Hex(1)HexNAc(1)NeuGc(2)
    composition.insert("[UNIMOD:1620]".to_string(), HashMap::from([("H", 59), ("C", 38), ("N", 3), ("O", 27)])); // Hex(1)HexNAc(1)NeuAc(2)Ac(1)
    composition.insert("[UNIMOD:1621]".to_string(), HashMap::from([("H", 61), ("C", 38), ("N", 1), ("O", 29)])); // dHex(2)Hex(2)HexA(1)HexNAc(1)
    composition.insert("[UNIMOD:1622]".to_string(), HashMap::from([("H", 59), ("C", 36), ("N", 3), ("O", 27), ("S", 1)])); // dHex(1)Hex(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1623]".to_string(), HashMap::from([("H", 53), ("C", 34), ("N", 1), ("O", 31), ("S", 1)])); // Hex(2)HexA(1)NeuAc(1)Pent(1)Sulf(1)
    composition.insert("[UNIMOD:1624]".to_string(), HashMap::from([("H", 63), ("C", 39), ("N", 3), ("O", 27)])); // dHex(1)Hex(1)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1625]".to_string(), HashMap::from([("H", 61), ("C", 38), ("N", 1), ("O", 30)])); // dHex(1)Hex(3)HexA(1)HexNAc(1)
    composition.insert("[UNIMOD:1626]".to_string(), HashMap::from([("H", 59), ("C", 36), ("N", 3), ("O", 28), ("S", 1)])); // Hex(2)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1627]".to_string(), HashMap::from([("H", 63), ("C", 38), ("N", 1), ("O", 30)])); // Hex(5)HexNAc(1)
    composition.insert("[UNIMOD:1628]".to_string(), HashMap::from([("H", 65), ("C", 40), ("N", 5), ("O", 25)])); // HexNAc(5)
    composition.insert("[UNIMOD:1630]".to_string(), HashMap::from([("H", 61), ("C", 40), ("N", 3), ("O", 28)])); // Hex(1)HexNAc(1)NeuAc(2)Ac(2)
    composition.insert("[UNIMOD:1631]".to_string(), HashMap::from([("H", 63), ("C", 39), ("N", 3), ("O", 29)])); // Hex(2)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1632]".to_string(), HashMap::from([("H", 53), ("C", 30), ("O", 34), ("P", 3)])); // Hex(5)Phos(3)
    composition.insert("[UNIMOD:1633]".to_string(), HashMap::from([("H", 61), ("C", 36), ("O", 33), ("P", 1)])); // Hex(6)Phos(1)
    composition.insert("[UNIMOD:1634]".to_string(), HashMap::from([("H", 64), ("C", 40), ("N", 2), ("O", 30)])); // dHex(1)Hex(2)HexA(1)HexNAc(2)
    composition.insert("[UNIMOD:1635]".to_string(), HashMap::from([("H", 63), ("C", 38), ("N", 1), ("O", 31), ("S", 1)])); // dHex(2)Hex(3)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1636]".to_string(), HashMap::from([("H", 66), ("C", 41), ("N", 4), ("O", 28)])); // Hex(1)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1637]".to_string(), HashMap::from([("H", 69), ("C", 42), ("N", 3), ("O", 28)])); // dHex(2)Hex(1)HexNAc(3)
    composition.insert("[UNIMOD:1638]".to_string(), HashMap::from([("H", 66), ("C", 41), ("N", 4), ("O", 29)])); // Hex(1)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1639]".to_string(), HashMap::from([("H", 63), ("C", 39), ("N", 3), ("O", 30), ("S", 1)])); // dHex(1)Hex(1)HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1640]".to_string(), HashMap::from([("H", 61), ("C", 38), ("N", 1), ("O", 33), ("S", 1)])); // dHex(1)Hex(3)HexA(1)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1641]".to_string(), HashMap::from([("H", 67), ("C", 42), ("N", 3), ("O", 30)])); // dHex(1)Hex(1)HexA(1)HexNAc(3)
    composition.insert("[UNIMOD:1642]".to_string(), HashMap::from([("H", 63), ("C", 39), ("N", 3), ("O", 31), ("S", 1)])); // Hex(2)HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1643]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 31), ("S", 1)])); // dHex(2)Hex(2)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1644]".to_string(), HashMap::from([("H", 70), ("C", 43), ("N", 2), ("O", 31)])); // dHex(2)Hex(1)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1645]".to_string(), HashMap::from([("H", 72), ("C", 44), ("N", 4), ("O", 29)])); // dHex(1)Hex(1)HexNAc(4)
    composition.insert("[UNIMOD:1646]".to_string(), HashMap::from([("H", 72), ("C", 44), ("N", 4), ("O", 30)])); // Hex(2)HexNAc(4)
    composition.insert("[UNIMOD:1647]".to_string(), HashMap::from([("H", 67), ("C", 42), ("N", 3), ("O", 33)])); // Hex(2)HexNAc(1)NeuGc(2)
    composition.insert("[UNIMOD:1648]".to_string(), HashMap::from([("H", 73), ("C", 44), ("N", 1), ("O", 33)])); // dHex(2)Hex(4)HexNAc(1)
    composition.insert("[UNIMOD:1649]".to_string(), HashMap::from([("H", 70), ("C", 44), ("N", 4), ("O", 31)])); // Hex(1)HexNAc(2)NeuAc(2)
    composition.insert("[UNIMOD:1650]".to_string(), HashMap::from([("H", 73), ("C", 45), ("N", 3), ("O", 31)])); // dHex(2)Hex(1)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1651]".to_string(), HashMap::from([("H", 69), ("C", 42), ("N", 3), ("O", 32), ("S", 1)])); // dHex(1)Hex(2)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1652]".to_string(), HashMap::from([("H", 75), ("C", 46), ("N", 5), ("O", 29)])); // dHex(1)HexNAc(5)
    composition.insert("[UNIMOD:1653]".to_string(), HashMap::from([("H", 73), ("C", 45), ("N", 3), ("O", 32)])); // dHex(2)Hex(1)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1654]".to_string(), HashMap::from([("H", 76), ("C", 46), ("N", 2), ("O", 32)])); // dHex(3)Hex(2)HexNAc(2)
    composition.insert("[UNIMOD:1655]".to_string(), HashMap::from([("H", 69), ("C", 42), ("N", 3), ("O", 33), ("S", 1)])); // Hex(3)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1656]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 34), ("S", 2)])); // dHex(2)Hex(2)HexNAc(2)Sulf(2)
    composition.insert("[UNIMOD:1657]".to_string(), HashMap::from([("H", 73), ("C", 45), ("N", 3), ("O", 33)])); // dHex(1)Hex(2)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1658]".to_string(), HashMap::from([("H", 76), ("C", 47), ("N", 4), ("O", 32)])); // dHex(1)Hex(1)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1659]".to_string(), HashMap::from([("H", 63), ("C", 36), ("O", 39), ("P", 3)])); // Hex(6)Phos(3)
    composition.insert("[UNIMOD:1660]".to_string(), HashMap::from([("H", 74), ("C", 46), ("N", 2), ("O", 35)])); // dHex(1)Hex(3)HexA(1)HexNAc(2)
    composition.insert("[UNIMOD:1661]".to_string(), HashMap::from([("H", 76), ("C", 47), ("N", 4), ("O", 33)])); // dHex(1)Hex(1)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1662]".to_string(), HashMap::from([("H", 70), ("C", 44), ("N", 4), ("O", 34), ("S", 1)])); // Hex(1)HexNAc(2)NeuAc(2)Sulf(1)
    composition.insert("[UNIMOD:1663]".to_string(), HashMap::from([("H", 71), ("C", 44), ("N", 1), ("O", 37), ("S", 1)])); // dHex(2)Hex(3)HexA(1)HexNAc(1)Sulf(1)
    composition.insert("[UNIMOD:1664]".to_string(), HashMap::from([("H", 74), ("C", 47), ("N", 4), ("O", 34)])); // Hex(1)HexNAc(1)NeuAc(3)
    composition.insert("[UNIMOD:1665]".to_string(), HashMap::from([("H", 76), ("C", 47), ("N", 4), ("O", 34)])); // Hex(2)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1666]".to_string(), HashMap::from([("H", 73), ("C", 45), ("N", 3), ("O", 35), ("S", 1)])); // dHex(1)Hex(2)HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1667]".to_string(), HashMap::from([("H", 80), ("C", 49), ("N", 2), ("O", 35)])); // dHex(3)Hex(1)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1668]".to_string(), HashMap::from([("H", 76), ("C", 46), ("N", 2), ("O", 36), ("S", 1)])); // dHex(2)Hex(3)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1669]".to_string(), HashMap::from([("H", 80), ("C", 49), ("N", 2), ("O", 36)])); // dHex(2)Hex(2)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1670]".to_string(), HashMap::from([("H", 74), ("C", 46), ("N", 2), ("O", 37), ("S", 1)])); // dHex(2)Hex(2)HexA(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1671]".to_string(), HashMap::from([("H", 82), ("C", 50), ("N", 4), ("O", 34)])); // dHex(1)Hex(2)HexNAc(4)
    composition.insert("[UNIMOD:1672]".to_string(), HashMap::from([("H", 74), ("C", 47), ("N", 4), ("O", 37)])); // Hex(1)HexNAc(1)NeuGc(3)
    composition.insert("[UNIMOD:1673]".to_string(), HashMap::from([("H", 76), ("C", 47), ("N", 4), ("O", 35), ("S", 1)])); // dHex(1)Hex(1)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1674]".to_string(), HashMap::from([("H", 74), ("C", 46), ("N", 2), ("O", 38), ("S", 1)])); // dHex(1)Hex(3)HexA(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1675]".to_string(), HashMap::from([("H", 80), ("C", 50), ("N", 4), ("O", 35)])); // dHex(1)Hex(1)HexNAc(2)NeuAc(2)
    composition.insert("[UNIMOD:1676]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 35)])); // dHex(3)HexNAc(3)Kdn(1)
    composition.insert("[UNIMOD:1678]".to_string(), HashMap::from([("H", 76), ("C", 47), ("N", 4), ("O", 36), ("S", 1)])); // Hex(2)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1679]".to_string(), HashMap::from([("H", 79), ("C", 48), ("N", 3), ("O", 36), ("S", 1)])); // dHex(2)Hex(2)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1680]".to_string(), HashMap::from([("H", 85), ("C", 52), ("N", 5), ("O", 33)])); // dHex(2)HexNAc(5)
    composition.insert("[UNIMOD:1681]".to_string(), HashMap::from([("H", 80), ("C", 50), ("N", 4), ("O", 36)])); // Hex(2)HexNAc(2)NeuAc(2)
    composition.insert("[UNIMOD:1682]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 36)])); // dHex(2)Hex(2)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1683]".to_string(), HashMap::from([("H", 79), ("C", 48), ("N", 3), ("O", 37), ("S", 1)])); // dHex(1)Hex(3)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1684]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 37)])); // dHex(2)Hex(2)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1685]".to_string(), HashMap::from([("H", 85), ("C", 52), ("N", 5), ("O", 35)])); // Hex(2)HexNAc(5)
    composition.insert("[UNIMOD:1686]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 38)])); // dHex(1)Hex(3)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1687]".to_string(), HashMap::from([("H", 83), ("C", 52), ("N", 5), ("O", 36)])); // Hex(1)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1688]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 37)])); // dHex(1)Hex(2)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1689]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 37)])); // dHex(3)Hex(2)HexNAc(3)
    composition.insert("[UNIMOD:1690]".to_string(), HashMap::from([("H", 73), ("C", 42), ("O", 44), ("P", 3)])); // Hex(7)Phos(3)
    composition.insert("[UNIMOD:1691]".to_string(), HashMap::from([("H", 84), ("C", 52), ("N", 2), ("O", 40)])); // dHex(1)Hex(4)HexA(1)HexNAc(2)
    composition.insert("[UNIMOD:1692]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 38)])); // Hex(3)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1693]".to_string(), HashMap::from([("H", 82), ("C", 52), ("N", 2), ("O", 41)])); // dHex(1)Hex(3)HexA(2)HexNAc(2)
    composition.insert("[UNIMOD:1694]".to_string(), HashMap::from([("H", 80), ("C", 50), ("N", 4), ("O", 39), ("S", 1)])); // Hex(2)HexNAc(2)NeuAc(2)Sulf(1)
    composition.insert("[UNIMOD:1695]".to_string(), HashMap::from([("H", 83), ("C", 51), ("N", 3), ("O", 39), ("S", 1)])); // dHex(2)Hex(2)HexNAc(2)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1696]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 39)])); // Hex(3)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1697]".to_string(), HashMap::from([("H", 90), ("C", 55), ("N", 2), ("O", 39)])); // dHex(4)Hex(1)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1698]".to_string(), HashMap::from([("H", 90), ("C", 55), ("N", 2), ("O", 40)])); // dHex(3)Hex(2)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1699]".to_string(), HashMap::from([("H", 84), ("C", 52), ("N", 2), ("O", 41), ("S", 1)])); // dHex(3)Hex(2)HexA(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1700]".to_string(), HashMap::from([("H", 89), ("C", 55), ("N", 5), ("O", 38)])); // Hex(2)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1701]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 38)])); // dHex(2)Hex(2)HexNAc(4)
    composition.insert("[UNIMOD:1702]".to_string(), HashMap::from([("H", 84), ("C", 52), ("N", 2), ("O", 42), ("S", 1)])); // dHex(2)Hex(3)HexA(1)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1703]".to_string(), HashMap::from([("H", 93), ("C", 57), ("N", 3), ("O", 39)])); // dHex(4)HexNAc(3)Kdn(1)
    composition.insert("[UNIMOD:1705]".to_string(), HashMap::from([("H", 84), ("C", 53), ("N", 4), ("O", 42)])); // Hex(2)HexNAc(1)NeuGc(3)
    composition.insert("[UNIMOD:1706]".to_string(), HashMap::from([("H", 91), ("C", 56), ("N", 1), ("O", 42)])); // dHex(4)Hex(1)HexNAc(1)Kdn(2)
    composition.insert("[UNIMOD:1707]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 40), ("S", 1)])); // dHex(1)Hex(2)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1708]".to_string(), HashMap::from([("H", 90), ("C", 56), ("N", 4), ("O", 40)])); // dHex(1)Hex(2)HexNAc(2)NeuAc(2)
    composition.insert("[UNIMOD:1709]".to_string(), HashMap::from([("H", 93), ("C", 57), ("N", 3), ("O", 40)])); // dHex(3)Hex(1)HexNAc(3)Kdn(1)
    composition.insert("[UNIMOD:1711]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 41), ("S", 1)])); // Hex(3)HexNAc(3)NeuAc(1)Sulf(1)
    composition.insert("[UNIMOD:1712]".to_string(), HashMap::from([("H", 90), ("C", 56), ("N", 4), ("O", 41)])); // Hex(3)HexNAc(2)NeuAc(2)
    composition.insert("[UNIMOD:1713]".to_string(), HashMap::from([("H", 86), ("C", 53), ("N", 4), ("O", 42), ("S", 1)])); // Hex(3)HexNAc(3)NeuGc(1)Sulf(1)
    composition.insert("[UNIMOD:1714]".to_string(), HashMap::from([("H", 90), ("C", 56), ("N", 4), ("O", 42)])); // dHex(1)Hex(2)HexNAc(2)NeuGc(2)
    composition.insert("[UNIMOD:1715]".to_string(), HashMap::from([("H", 93), ("C", 57), ("N", 3), ("O", 42)])); // dHex(2)Hex(3)HexNAc(2)NeuGc(1)
    composition.insert("[UNIMOD:1716]".to_string(), HashMap::from([("H", 87), ("C", 54), ("N", 3), ("O", 43), ("S", 1)])); // dHex(1)Hex(3)HexA(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1717]".to_string(), HashMap::from([("H", 93), ("C", 58), ("N", 5), ("O", 41)])); // Hex(2)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1718]".to_string(), HashMap::from([("H", 96), ("C", 59), ("N", 4), ("O", 41)])); // dHex(2)Hex(2)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1719]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 41)])); // dHex(4)Hex(2)HexNAc(3)
    composition.insert("[UNIMOD:1720]".to_string(), HashMap::from([("H", 93), ("C", 58), ("N", 5), ("O", 42)])); // Hex(2)HexNAc(3)NeuAc(1)NeuGc(1)
    composition.insert("[UNIMOD:1721]".to_string(), HashMap::from([("H", 96), ("C", 59), ("N", 4), ("O", 42)])); // dHex(2)Hex(2)HexNAc(3)NeuGc(1)
    composition.insert("[UNIMOD:1722]".to_string(), HashMap::from([("H", 99), ("C", 60), ("N", 3), ("O", 42)])); // dHex(3)Hex(3)HexNAc(3)
    composition.insert("[UNIMOD:1723]".to_string(), HashMap::from([("H", 83), ("C", 48), ("O", 49), ("P", 3)])); // Hex(8)Phos(3)
    composition.insert("[UNIMOD:1724]".to_string(), HashMap::from([("H", 90), ("C", 56), ("N", 4), ("O", 43), ("S", 1)])); // dHex(1)Hex(2)HexNAc(2)NeuAc(2)Sulf(1)
    composition.insert("[UNIMOD:1725]".to_string(), HashMap::from([("H", 93), ("C", 58), ("N", 5), ("O", 43)])); // Hex(2)HexNAc(3)NeuGc(2)
    composition.insert("[UNIMOD:1726]".to_string(), HashMap::from([("H", 100), ("C", 61), ("N", 2), ("O", 44)])); // dHex(4)Hex(2)HexNAc(2)Kdn(1)
    composition.insert("[UNIMOD:1727]".to_string(), HashMap::from([("H", 99), ("C", 61), ("N", 5), ("O", 42)])); // dHex(1)Hex(2)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1728]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 42)])); // dHex(3)Hex(2)HexNAc(4)
    composition.insert("[UNIMOD:1729]".to_string(), HashMap::from([("H", 91), ("C", 58), ("N", 5), ("O", 46)])); // Hex(1)HexNAc(1)NeuGc(4)
    composition.insert("[UNIMOD:1730]".to_string(), HashMap::from([("H", 103), ("C", 63), ("N", 3), ("O", 44)])); // dHex(4)Hex(1)HexNAc(3)Kdn(1)
    composition.insert("[UNIMOD:1732]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 46), ("S", 2)])); // Hex(4)HexNAc(4)Sulf(2)
    composition.insert("[UNIMOD:1733]".to_string(), HashMap::from([("H", 103), ("C", 63), ("N", 3), ("O", 45)])); // dHex(3)Hex(2)HexNAc(3)Kdn(1)
    composition.insert("[UNIMOD:1735]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 5), ("O", 43)])); // dHex(2)Hex(2)HexNAc(5)
    composition.insert("[UNIMOD:1736]".to_string(), HashMap::from([("H", 97), ("C", 60), ("N", 3), ("O", 47), ("S", 1)])); // dHex(2)Hex(3)HexA(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1737]".to_string(), HashMap::from([("H", 97), ("C", 60), ("N", 3), ("O", 48), ("S", 1)])); // dHex(1)Hex(4)HexA(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1738]".to_string(), HashMap::from([("H", 103), ("C", 64), ("N", 5), ("O", 46)])); // Hex(3)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1739]".to_string(), HashMap::from([("H", 106), ("C", 65), ("N", 4), ("O", 46)])); // dHex(2)Hex(3)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1740]".to_string(), HashMap::from([("H", 109), ("C", 66), ("N", 3), ("O", 46)])); // dHex(4)Hex(3)HexNAc(3)
    composition.insert("[UNIMOD:1742]".to_string(), HashMap::from([("H", 93), ("C", 54), ("O", 54), ("P", 3)])); // Hex(9)Phos(3)
    composition.insert("[UNIMOD:1743]".to_string(), HashMap::from([("H", 111), ("C", 68), ("N", 7), ("O", 43)])); // dHex(2)HexNAc(7)
    composition.insert("[UNIMOD:1744]".to_string(), HashMap::from([("H", 101), ("C", 64), ("N", 5), ("O", 51)])); // Hex(2)HexNAc(1)NeuGc(4)
    composition.insert("[UNIMOD:1745]".to_string(), HashMap::from([("H", 103), ("C", 64), ("N", 5), ("O", 49), ("S", 1)])); // Hex(3)HexNAc(3)NeuAc(2)Sulf(1)
    composition.insert("[UNIMOD:1746]".to_string(), HashMap::from([("H", 115), ("C", 70), ("N", 5), ("O", 48)])); // dHex(2)Hex(3)HexNAc(5)
    composition.insert("[UNIMOD:1747]".to_string(), HashMap::from([("H", 107), ("C", 67), ("N", 5), ("O", 51)])); // dHex(1)Hex(2)HexNAc(2)NeuGc(3)
    composition.insert("[UNIMOD:1748]".to_string(), HashMap::from([("H", 107), ("C", 66), ("N", 3), ("O", 52), ("S", 1)])); // dHex(2)Hex(4)HexA(1)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1749]".to_string(), HashMap::from([("H", 110), ("C", 69), ("N", 6), ("O", 49)])); // Hex(2)HexNAc(3)NeuAc(3)
    composition.insert("[UNIMOD:1750]".to_string(), HashMap::from([("H", 113), ("C", 70), ("N", 5), ("O", 50)])); // dHex(1)Hex(3)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1751]".to_string(), HashMap::from([("H", 116), ("C", 71), ("N", 4), ("O", 50)])); // dHex(3)Hex(3)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1752]".to_string(), HashMap::from([("H", 110), ("C", 69), ("N", 6), ("O", 52)])); // Hex(2)HexNAc(3)NeuGc(3)
    composition.insert("[UNIMOD:1753]".to_string(), HashMap::from([("H", 103), ("C", 60), ("O", 59), ("P", 3)])); // Hex(10)Phos(3)
    composition.insert("[UNIMOD:1754]".to_string(), HashMap::from([("H", 116), ("C", 72), ("N", 6), ("O", 50)])); // dHex(1)Hex(2)HexNAc(4)NeuAc(2)
    composition.insert("[UNIMOD:1755]".to_string(), HashMap::from([("H", 108), ("C", 69), ("N", 6), ("O", 55)])); // Hex(1)HexNAc(1)NeuGc(5)
    composition.insert("[UNIMOD:1756]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 54), ("S", 2)])); // Hex(4)HexNAc(4)NeuAc(1)Sulf(2)
    composition.insert("[UNIMOD:1757]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 55), ("S", 2)])); // Hex(4)HexNAc(4)NeuGc(1)Sulf(2)
    composition.insert("[UNIMOD:1758]".to_string(), HashMap::from([("H", 123), ("C", 76), ("N", 5), ("O", 54)])); // dHex(2)Hex(3)HexNAc(3)NeuAc(2)
    composition.insert("[UNIMOD:1759]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 57), ("S", 3)])); // Hex(4)HexNAc(4)NeuAc(1)Sulf(3)
    composition.insert("[UNIMOD:1760]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 28)])); // dHex(2)Hex(2)HexNAc(2)
    composition.insert("[UNIMOD:1761]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 29)])); // dHex(1)Hex(3)HexNAc(2)
    composition.insert("[UNIMOD:1762]".to_string(), HashMap::from([("H", 69), ("C", 42), ("N", 3), ("O", 29)])); // dHex(1)Hex(2)HexNAc(3)
    composition.insert("[UNIMOD:1763]".to_string(), HashMap::from([("H", 69), ("C", 42), ("N", 3), ("O", 30)])); // Hex(3)HexNAc(3)
    composition.insert("[UNIMOD:1764]".to_string(), HashMap::from([("H", 66), ("C", 40), ("N", 2), ("O", 32), ("S", 1)])); // dHex(1)Hex(3)HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1765]".to_string(), HashMap::from([("H", 76), ("C", 46), ("N", 2), ("O", 33)])); // dHex(2)Hex(3)HexNAc(2)
    composition.insert("[UNIMOD:1766]".to_string(), HashMap::from([("H", 76), ("C", 46), ("N", 2), ("O", 34)])); // dHex(1)Hex(4)HexNAc(2)
    composition.insert("[UNIMOD:1767]".to_string(), HashMap::from([("H", 79), ("C", 48), ("N", 3), ("O", 33)])); // dHex(2)Hex(2)HexNAc(3)
    composition.insert("[UNIMOD:1768]".to_string(), HashMap::from([("H", 79), ("C", 48), ("N", 3), ("O", 34)])); // dHex(1)Hex(3)HexNAc(3)
    composition.insert("[UNIMOD:1769]".to_string(), HashMap::from([("H", 79), ("C", 48), ("N", 3), ("O", 35)])); // Hex(4)HexNAc(3)
    composition.insert("[UNIMOD:1770]".to_string(), HashMap::from([("H", 86), ("C", 52), ("N", 2), ("O", 38)])); // dHex(2)Hex(4)HexNAc(2)
    composition.insert("[UNIMOD:1771]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 38)])); // dHex(2)Hex(3)HexNAc(3)
    composition.insert("[UNIMOD:1772]".to_string(), HashMap::from([("H", 95), ("C", 58), ("N", 5), ("O", 40)])); // Hex(3)HexNAc(5)
    composition.insert("[UNIMOD:1773]".to_string(), HashMap::from([("H", 96), ("C", 59), ("N", 4), ("O", 43)])); // Hex(4)HexNAc(3)NeuAc(1)
    composition.insert("[UNIMOD:1774]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 43)])); // dHex(2)Hex(3)HexNAc(4)
    composition.insert("[UNIMOD:1775]".to_string(), HashMap::from([("H", 105), ("C", 64), ("N", 5), ("O", 44)])); // dHex(1)Hex(3)HexNAc(5)
    composition.insert("[UNIMOD:1776]".to_string(), HashMap::from([("H", 108), ("C", 66), ("N", 6), ("O", 45)])); // Hex(3)HexNAc(6)
    composition.insert("[UNIMOD:1777]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 48)])); // Hex(4)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1778]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 48)])); // dHex(2)Hex(4)HexNAc(4)
    composition.insert("[UNIMOD:1779]".to_string(), HashMap::from([("H", 112), ("C", 68), ("N", 4), ("O", 50)])); // Hex(6)HexNAc(4)
    composition.insert("[UNIMOD:1780]".to_string(), HashMap::from([("H", 115), ("C", 70), ("N", 5), ("O", 50)])); // Hex(5)HexNAc(5)
    composition.insert("[UNIMOD:1781]".to_string(), HashMap::from([("H", 118), ("C", 72), ("N", 6), ("O", 49)])); // dHex(1)Hex(3)HexNAc(6)
    composition.insert("[UNIMOD:1782]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 52)])); // dHex(1)Hex(4)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1783]".to_string(), HashMap::from([("H", 122), ("C", 74), ("N", 4), ("O", 52)])); // dHex(3)Hex(4)HexNAc(4)
    composition.insert("[UNIMOD:1784]".to_string(), HashMap::from([("H", 122), ("C", 75), ("N", 6), ("O", 52)])); // dHex(1)Hex(3)HexNAc(5)NeuAc(1)
    composition.insert("[UNIMOD:1785]".to_string(), HashMap::from([("H", 125), ("C", 76), ("N", 5), ("O", 53)])); // dHex(2)Hex(4)HexNAc(5)
    composition.insert("[UNIMOD:1786]".to_string(), HashMap::from([("H", 42), ("C", 27), ("N", 2), ("O", 19)])); // Hex(1)HexNAc(1)NeuAc(1)Ac(1)
    composition.insert("[UNIMOD:1787]".to_string(), HashMap::from([("C", -2), ("13C", 2), ("N", -2), ("15N", 2)])); // Label:13C(2)15N(2)
    composition.insert("[UNIMOD:1789]".to_string(), HashMap::from([("H", 13), ("C", 8), ("N", 1), ("O", 2)])); // Xlink:DSS[155]
    composition.insert("[UNIMOD:1799]".to_string(), HashMap::from([("H", 31), ("C", 19), ("N", 7), ("O", 7)])); // NQIGG
    composition.insert("[UNIMOD:1800]".to_string(), HashMap::from([("H", 6), ("C", 7), ("O", 2)])); // Carboxyethylpyrrole
    composition.insert("[UNIMOD:1801]".to_string(), HashMap::from([("H", 19), ("C", 29), ("N", 1), ("O", 7)])); // Fluorescein-tyramine
    composition.insert("[UNIMOD:1824]".to_string(), HashMap::from([("H", 6), ("C", 4), ("O", 2)])); // GEE
    composition.insert("[UNIMOD:1825]".to_string(), HashMap::from([("H", 13), ("C", 9), ("N", 2), ("O", 9), ("P", 1)])); // RNPXL
    composition.insert("[UNIMOD:1826]".to_string(), HashMap::from([("C", 1), ("O", -1)])); // Glu->pyro-Glu+Methyl
    composition.insert("[UNIMOD:1827]".to_string(), HashMap::from([("H", -2), ("2H", 2), ("13C", 1), ("O", -1)])); // Glu->pyro-Glu+Methyl:2H(2)13C(1)
    composition.insert("[UNIMOD:1828]".to_string(), HashMap::from([("H", 31), ("C", 17), ("N", 7), ("O", 4)])); // LRGG+methyl
    composition.insert("[UNIMOD:1829]".to_string(), HashMap::from([("H", 33), ("C", 18), ("N", 7), ("O", 4)])); // LRGG+dimethyl
    composition.insert("[UNIMOD:1830]".to_string(), HashMap::from([("H", 23), ("C", 18), ("N", 3), ("O", 3), ("S", 1)])); // Biotin-tyramide
    composition.insert("[UNIMOD:1831]".to_string(), HashMap::from([("H", 10), ("C", 4), ("N", 1), ("O", 2)])); // Tris
    composition.insert("[UNIMOD:1832]".to_string(), HashMap::from([("H", 16), ("C", 18), ("N", 2), ("O", 8), ("S", 2)])); // IASD
    composition.insert("[UNIMOD:1833]".to_string(), HashMap::from([("H", 24), ("C", 15), ("O", 1)])); // NP40
    composition.insert("[UNIMOD:1834]".to_string(), HashMap::from([("H", 21), ("C", 12)])); // Tween20
    composition.insert("[UNIMOD:1835]".to_string(), HashMap::from([("H", 31), ("C", 18), ("O", 1)])); // Tween80
    composition.insert("[UNIMOD:1836]".to_string(), HashMap::from([("H", 20), ("C", 14)])); // Triton
    composition.insert("[UNIMOD:1837]".to_string(), HashMap::from([("H", 24), ("C", 12)])); // Brij35
    composition.insert("[UNIMOD:1838]".to_string(), HashMap::from([("H", 32), ("C", 16)])); // Brij58
    composition.insert("[UNIMOD:1839]".to_string(), HashMap::from([("H", 30), ("C", 25), ("N", 2), ("O", 6)])); // betaFNA
    composition.insert("[UNIMOD:1840]".to_string(), HashMap::from([("H", 132), ("C", 80), ("N", 4), ("O", 59)])); // dHex(1)Hex(7)HexNAc(4)
    composition.insert("[UNIMOD:1841]".to_string(), HashMap::from([("H", 23), ("C", 15), ("N", 3), ("O", 3), ("S", 3)])); // Biotin:Thermo-21328
    composition.insert("[UNIMOD:1843]".to_string(), HashMap::from([("H", 12), ("C", 9), ("N", 3), ("O", 7), ("P", 1)])); // PhosphoCytidine
    composition.insert("[UNIMOD:1845]".to_string(), HashMap::from([("H", -1), ("N", 3)])); // AzidoF
    composition.insert("[UNIMOD:1846]".to_string(), HashMap::from([("H", 9), ("C", 4), ("N", 1)])); // Dimethylaminoethyl
    composition.insert("[UNIMOD:1848]".to_string(), HashMap::from([("H", 6), ("C", 5), ("O", 3)])); // Gluratylation
    composition.insert("[UNIMOD:1849]".to_string(), HashMap::from([("H", 6), ("C", 4), ("O", 2)])); // hydroxyisobutyryl
    composition.insert("[UNIMOD:1868]".to_string(), HashMap::from([("H", 5), ("C", 2), ("O", 1), ("P", 1), ("S", 1)])); // MeMePhosphorothioate
    composition.insert("[UNIMOD:1870]".to_string(), HashMap::from([("H", -3), ("Fe", 1)])); // Cation:Fe[III]
    composition.insert("[UNIMOD:1871]".to_string(), HashMap::from([("H", 8), ("C", 4), ("O", 2), ("S", 2)])); // DTT
    composition.insert("[UNIMOD:1872]".to_string(), HashMap::from([("H", 13), ("C", 11), ("O", 1)])); // DYn-2
    composition.insert("[UNIMOD:1873]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 1)])); // MesitylOxide
    composition.insert("[UNIMOD:1875]".to_string(), HashMap::from([("H", 2), ("C", 1), ("O", 1)])); // methylol
    composition.insert("[UNIMOD:1877]".to_string(), HashMap::from([("H", 21), ("C", 12), ("N", 1), ("O", 5)])); // Xlink:DSS[259]
    composition.insert("[UNIMOD:1878]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 4), ("S", 1)])); // Xlink:DSSO[176]
    composition.insert("[UNIMOD:1879]".to_string(), HashMap::from([("H", 9), ("C", 6), ("N", 1), ("O", 3), ("S", 1)])); // Xlink:DSSO[175]
    composition.insert("[UNIMOD:1880]".to_string(), HashMap::from([("H", 17), ("C", 10), ("N", 1), ("O", 6), ("S", 1)])); // Xlink:DSSO[279]
    composition.insert("[UNIMOD:1881]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 1)])); // Xlink:DSSO[54]
    composition.insert("[UNIMOD:1882]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 1), ("S", 1)])); // Xlink:DSSO[86]
    composition.insert("[UNIMOD:1883]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 2), ("S", 1)])); // Xlink:DSSO[104]
    composition.insert("[UNIMOD:1885]".to_string(), HashMap::from([("H", 5), ("C", 5), ("N", 1), ("O", 2)])); // Xlink:BuUrBu[111]
    composition.insert("[UNIMOD:1886]".to_string(), HashMap::from([("H", 7), ("C", 4), ("N", 1), ("O", 1)])); // Xlink:BuUrBu[85]
    composition.insert("[UNIMOD:1887]".to_string(), HashMap::from([("H", 15), ("C", 9), ("N", 3), ("O", 3)])); // Xlink:BuUrBu[213]
    composition.insert("[UNIMOD:1888]".to_string(), HashMap::from([("H", 14), ("C", 9), ("N", 2), ("O", 4)])); // Xlink:BuUrBu[214]
    composition.insert("[UNIMOD:1889]".to_string(), HashMap::from([("H", 23), ("C", 13), ("N", 3), ("O", 6)])); // Xlink:BuUrBu[317]
    composition.insert("[UNIMOD:1896]".to_string(), HashMap::from([("H", 6), ("C", 6), ("O", 3), ("S", 1)])); // Xlink:DSSO[158]
    composition.insert("[UNIMOD:1897]".to_string(), HashMap::from([("H", 10), ("C", 10), ("O", 6)])); // Xlink:EGS[226]
    composition.insert("[UNIMOD:1898]".to_string(), HashMap::from([("H", 10), ("C", 8), ("O", 2)])); // Xlink:DSS[138]
    composition.insert("[UNIMOD:1899]".to_string(), HashMap::from([("H", 12), ("C", 9), ("N", 2), ("O", 3)])); // Xlink:BuUrBu[196]
    composition.insert("[UNIMOD:1900]".to_string(), HashMap::from([("H", 8), ("C", 6), ("N", 2), ("S", 2)])); // Xlink:DTBP[172]
    composition.insert("[UNIMOD:1901]".to_string(), HashMap::from([("H", 2), ("C", 4), ("O", 4)])); // Xlink:DST[114]
    composition.insert("[UNIMOD:1902]".to_string(), HashMap::from([("H", 6), ("C", 6), ("O", 2), ("S", 2)])); // Xlink:DTSSP[174]
    composition.insert("[UNIMOD:1903]".to_string(), HashMap::from([("H", 13), ("C", 12), ("N", 1), ("O", 3)])); // Xlink:SMCC[219]
    composition.insert("[UNIMOD:1905]".to_string(), HashMap::from([("H", 4), ("C", 5), ("O", 2)])); // Xlink:BS2G[96]
    composition.insert("[UNIMOD:1906]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 2)])); // Xlink:BS2G[113]
    composition.insert("[UNIMOD:1907]".to_string(), HashMap::from([("H", 6), ("C", 5), ("O", 3)])); // Xlink:BS2G[114]
    composition.insert("[UNIMOD:1908]".to_string(), HashMap::from([("H", 15), ("C", 9), ("N", 1), ("O", 5)])); // Xlink:BS2G[217]
    composition.insert("[UNIMOD:1910]".to_string(), HashMap::from([("H", -3), ("Al", 1)])); // Cation:Al[III]
    composition.insert("[UNIMOD:1911]".to_string(), HashMap::from([("H", 13), ("C", 7), ("N", 3)])); // Xlink:DMP[139]
    composition.insert("[UNIMOD:1912]".to_string(), HashMap::from([("H", 10), ("C", 7), ("N", 2)])); // Xlink:DMP[122]
    composition.insert("[UNIMOD:1913]".to_string(), HashMap::from([("H", -2), ("C", 2)])); // glyoxalAGE
    composition.insert("[UNIMOD:1914]".to_string(), HashMap::from([("H", -4), ("C", -1), ("O", 1), ("S", -1)])); // Met->AspSA
    composition.insert("[UNIMOD:1915]".to_string(), HashMap::from([("H", -2), ("C", -1), ("O", -1)])); // Decarboxylation
    composition.insert("[UNIMOD:1916]".to_string(), HashMap::from([("H", -2), ("C", -1), ("N", -2), ("O", 2)])); // Aspartylurea
    composition.insert("[UNIMOD:1917]".to_string(), HashMap::from([("H", -1), ("C", -1), ("N", -1), ("O", 2)])); // Formylasparagine
    composition.insert("[UNIMOD:1918]".to_string(), HashMap::from([("H", -2), ("O", 1)])); // Carbonyl
    composition.insert("[UNIMOD:1920]".to_string(), HashMap::from([("H", 10), ("C", 17), ("O", 6)])); // AFB1_Dialdehyde
    composition.insert("[UNIMOD:1922]".to_string(), HashMap::from([("H", 2), ("O", 1)])); // Pro->HAVA
    composition.insert("[UNIMOD:1923]".to_string(), HashMap::from([("H", -4), ("O", 2)])); // Delta:H(-4)O(2)
    composition.insert("[UNIMOD:1924]".to_string(), HashMap::from([("H", -4), ("O", 3)])); // Delta:H(-4)O(3)
    composition.insert("[UNIMOD:1925]".to_string(), HashMap::from([("O", 4)])); // Delta:O(4)
    composition.insert("[UNIMOD:1926]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 2)])); // Delta:H(4)C(3)O(2)
    composition.insert("[UNIMOD:1927]".to_string(), HashMap::from([("H", 4), ("C", 5), ("O", 1)])); // Delta:H(4)C(5)O(1)
    composition.insert("[UNIMOD:1928]".to_string(), HashMap::from([("H", 10), ("C", 8), ("O", 1)])); // Delta:H(10)C(8)O(1)
    composition.insert("[UNIMOD:1929]".to_string(), HashMap::from([("H", 6), ("C", 7), ("O", 4)])); // Delta:H(6)C(7)O(4)
    composition.insert("[UNIMOD:1930]".to_string(), HashMap::from([("H", 16), ("C", 10), ("O", 8)])); // Pent(2)
    composition.insert("[UNIMOD:1931]".to_string(), HashMap::from([("H", 21), ("C", 13), ("N", 1), ("O", 9)])); // Pent(1)HexNAc(1)
    composition.insert("[UNIMOD:1932]".to_string(), HashMap::from([("H", 20), ("C", 12), ("O", 13), ("S", 1)])); // Hex(2)Sulf(1)
    composition.insert("[UNIMOD:1933]".to_string(), HashMap::from([("H", 28), ("C", 17), ("O", 13)])); // Hex(1)Pent(2)Me(1)
    composition.insert("[UNIMOD:1934]".to_string(), HashMap::from([("H", 26), ("C", 16), ("N", 2), ("O", 13), ("S", 1)])); // HexNAc(2)Sulf(1)
    composition.insert("[UNIMOD:1935]".to_string(), HashMap::from([("H", 36), ("C", 22), ("O", 17)])); // Hex(1)Pent(3)Me(1)
    composition.insert("[UNIMOD:1936]".to_string(), HashMap::from([("H", 36), ("C", 22), ("O", 18)])); // Hex(2)Pent(2)
    composition.insert("[UNIMOD:1937]".to_string(), HashMap::from([("H", 38), ("C", 23), ("O", 18)])); // Hex(2)Pent(2)Me(1)
    composition.insert("[UNIMOD:1938]".to_string(), HashMap::from([("H", 48), ("C", 30), ("O", 26)])); // Hex(4)HexA(1)
    composition.insert("[UNIMOD:1939]".to_string(), HashMap::from([("H", 49), ("C", 31), ("N", 1), ("O", 25)])); // Hex(2)HexNAc(1)Pent(1)HexA(1)
    composition.insert("[UNIMOD:1940]".to_string(), HashMap::from([("H", 51), ("C", 32), ("N", 1), ("O", 26)])); // Hex(3)HexNAc(1)HexA(1)
    composition.insert("[UNIMOD:1941]".to_string(), HashMap::from([("H", 56), ("C", 34), ("N", 2), ("O", 26), ("S", 1)])); // Hex(1)HexNAc(2)dHex(2)Sulf(1)
    composition.insert("[UNIMOD:1942]".to_string(), HashMap::from([("H", 55), ("C", 36), ("N", 3), ("O", 27)])); // HexA(2)HexNAc(3)
    composition.insert("[UNIMOD:1943]".to_string(), HashMap::from([("H", 58), ("C", 36), ("O", 30)])); // dHex(1)Hex(4)HexA(1)
    composition.insert("[UNIMOD:1944]".to_string(), HashMap::from([("H", 58), ("C", 36), ("O", 31)])); // Hex(5)HexA(1)
    composition.insert("[UNIMOD:1945]".to_string(), HashMap::from([("H", 61), ("C", 38), ("N", 1), ("O", 31)])); // Hex(4)HexA(1)HexNAc(1)
    composition.insert("[UNIMOD:1946]".to_string(), HashMap::from([("H", 73), ("C", 44), ("N", 1), ("O", 32)])); // dHex(3)Hex(3)HexNAc(1)
    composition.insert("[UNIMOD:1947]".to_string(), HashMap::from([("H", 73), ("C", 44), ("N", 1), ("O", 35)])); // Hex(6)HexNAc(1)
    composition.insert("[UNIMOD:1948]".to_string(), HashMap::from([("H", 72), ("C", 44), ("N", 4), ("O", 32), ("S", 1)])); // Hex(1)HexNAc(4)dHex(1)Sulf(1)
    composition.insert("[UNIMOD:1949]".to_string(), HashMap::from([("H", 77), ("C", 48), ("N", 3), ("O", 35)])); // dHex(1)Hex(2)HexNAc(1)NeuAc(2)
    composition.insert("[UNIMOD:1950]".to_string(), HashMap::from([("H", 86), ("C", 52), ("N", 2), ("O", 37)])); // dHex(3)Hex(3)HexNAc(2)
    composition.insert("[UNIMOD:1951]".to_string(), HashMap::from([("H", 82), ("C", 50), ("N", 4), ("O", 36), ("S", 1)])); // dHex(2)Hex(1)HexNAc(4)Sulf(1)
    composition.insert("[UNIMOD:1952]".to_string(), HashMap::from([("H", 82), ("C", 50), ("N", 4), ("O", 40), ("S", 2)])); // dHex(1)Hex(2)HexNAc(4)Sulf(2)
    composition.insert("[UNIMOD:1953]".to_string(), HashMap::from([("H", 90), ("C", 54), ("O", 45)])); // Hex(9)
    composition.insert("[UNIMOD:1954]".to_string(), HashMap::from([("H", 89), ("C", 54), ("N", 3), ("O", 41), ("S", 1)])); // dHex(2)Hex(3)HexNAc(3)Sulf(1)
    composition.insert("[UNIMOD:1955]".to_string(), HashMap::from([("H", 98), ("C", 59), ("N", 2), ("O", 43)])); // dHex(2)Hex(5)HexNAc(2)Me(1)
    composition.insert("[UNIMOD:1956]".to_string(), HashMap::from([("H", 92), ("C", 56), ("N", 4), ("O", 44), ("S", 2)])); // dHex(2)Hex(2)HexNAc(4)Sulf(2)
    composition.insert("[UNIMOD:1957]".to_string(), HashMap::from([("H", 103), ("C", 62), ("N", 1), ("O", 50)])); // Hex(9)HexNAc(1)
    composition.insert("[UNIMOD:1958]".to_string(), HashMap::from([("H", 102), ("C", 62), ("N", 4), ("O", 48), ("S", 2)])); // dHex(3)Hex(2)HexNAc(4)Sulf(2)
    composition.insert("[UNIMOD:1959]".to_string(), HashMap::from([("H", 109), ("C", 67), ("N", 5), ("O", 49)])); // Hex(4)HexNAc(4)NeuGc(1)
    composition.insert("[UNIMOD:1960]".to_string(), HashMap::from([("H", 113), ("C", 69), ("N", 3), ("O", 49)])); // dHex(4)Hex(3)HexNAc(2)NeuAc(1)
    composition.insert("[UNIMOD:1961]".to_string(), HashMap::from([("H", 112), ("C", 69), ("N", 6), ("O", 48)])); // Hex(3)HexNAc(5)NeuAc(1)
    composition.insert("[UNIMOD:1962]".to_string(), HashMap::from([("H", 113), ("C", 68), ("N", 1), ("O", 55)])); // Hex(10)HexNAc(1)
    composition.insert("[UNIMOD:1963]".to_string(), HashMap::from([("H", 116), ("C", 70), ("N", 2), ("O", 54)])); // dHex(1)Hex(8)HexNAc(2)
    composition.insert("[UNIMOD:1964]".to_string(), HashMap::from([("H", 116), ("C", 72), ("N", 6), ("O", 51)])); // Hex(3)HexNAc(4)NeuAc(2)
    composition.insert("[UNIMOD:1965]".to_string(), HashMap::from([("H", 119), ("C", 73), ("N", 5), ("O", 51)])); // dHex(2)Hex(3)HexNAc(4)NeuAc(1)
    composition.insert("[UNIMOD:1966]".to_string(), HashMap::from([("H", 118), ("C", 72), ("N", 6), ("O", 51), ("S", 1)])); // dHex(2)Hex(2)HexNAc(6)Sulf(1)
    composition.insert("[UNIMOD:1967]".to_string(), HashMap::from([("H", 121), ("C", 75), ("N", 5), ("O", 54)])); // Hex(5)HexNAc(4)NeuAc(1)Ac(1)
    composition.insert("[UNIMOD:1968]".to_string(), HashMap::from([("H", 120), ("C", 75), ("N", 6), ("O", 54)])); // Hex(3)HexNAc(3)NeuAc(3)
    composition.insert("[UNIMOD:1969]".to_string(), HashMap::from([("H", 123), ("C", 77), ("N", 5), ("O", 55)])); // Hex(5)HexNAc(4)NeuAc(1)Ac(2)
    composition.insert("[UNIMOD:1970]".to_string(), HashMap::from([("H", 18), ("C", 8), ("O", 3)])); // Unknown:162
    composition.insert("[UNIMOD:1971]".to_string(), HashMap::from([("H", -7), ("O", 1), ("Fe", 3)])); // Unknown:177
    composition.insert("[UNIMOD:1972]".to_string(), HashMap::from([("H", 22), ("C", 13), ("O", 2)])); // Unknown:210
    composition.insert("[UNIMOD:1973]".to_string(), HashMap::from([("H", 16), ("C", 10), ("O", 5)])); // Unknown:216
    composition.insert("[UNIMOD:1974]".to_string(), HashMap::from([("H", 14), ("C", 9), ("O", 7)])); // Unknown:234
    composition.insert("[UNIMOD:1975]".to_string(), HashMap::from([("H", 28), ("C", 13), ("O", 4)])); // Unknown:248
    composition.insert("[UNIMOD:1976]".to_string(), HashMap::from([("H", 4), ("C", 10), ("N", 1), ("O", 5), ("S", 1)])); // Unknown:250
    composition.insert("[UNIMOD:1977]".to_string(), HashMap::from([("H", 8), ("C", 4), ("N", 5), ("O", 7), ("S", 2)])); // Unknown:302
    composition.insert("[UNIMOD:1978]".to_string(), HashMap::from([("H", 18), ("C", 12), ("O", 9)])); // Unknown:306
    composition.insert("[UNIMOD:1979]".to_string(), HashMap::from([("H", 24), ("C", 12), ("N", 2), ("O", 6), ("S", 4)])); // Unknown:420
    composition.insert("[UNIMOD:1986]".to_string(), HashMap::from([("H", 9), ("C", 4), ("O", 2), ("P", 1), ("S", 1)])); // Diethylphosphothione
    composition.insert("[UNIMOD:1987]".to_string(), HashMap::from([("H", 5), ("C", 2), ("O", 2), ("P", 1), ("S", 1)])); // Dimethylphosphothione
    composition.insert("[UNIMOD:1989]".to_string(), HashMap::from([("H", 3), ("C", 1), ("O", 2), ("P", 1), ("S", 1)])); // monomethylphosphothione
    composition.insert("[UNIMOD:1990]".to_string(), HashMap::from([("H", 22), ("C", 13), ("N", 4), ("O", 4), ("S", 1)])); // CIGG
    composition.insert("[UNIMOD:1991]".to_string(), HashMap::from([("H", 92), ("C", 61), ("N", 14), ("O", 15), ("S", 2)])); // GNLLFLACYCIGG
    composition.insert("[UNIMOD:1992]".to_string(), HashMap::from([("H", 9), ("C", 10), ("N", 1), ("O", 1)])); // serotonylation
    composition.insert("[UNIMOD:1993]".to_string(), HashMap::from([("H", 33), ("C", 20), ("13C", 9), ("O", 10), ("P", 1)])); // TMPP-Ac:13C(9)
    composition.insert("[UNIMOD:1999]".to_string(), HashMap::from([("C", 2), ("O", 2)])); // Xlink:DST[56]
    composition.insert("[UNIMOD:2001]".to_string(), HashMap::from([("H", 16), ("C", 15), ("N", 2), ("O", 6)])); // ZQG
    composition.insert("[UNIMOD:2006]".to_string(), HashMap::from([("H", 7), ("C", 4), ("O", 3), ("P", 1), ("Cl", 2)])); // Haloxon
    composition.insert("[UNIMOD:2007]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 1), ("O", 1), ("P", 1), ("S", 1)])); // Methamidophos-S
    composition.insert("[UNIMOD:2008]".to_string(), HashMap::from([("H", 4), ("C", 1), ("N", 1), ("O", 2), ("P", 1)])); // Methamidophos-O
    composition.insert("[UNIMOD:2014]".to_string(), HashMap::from([("H", -1), ("N", 1)])); // Nitrene
    composition.insert("[UNIMOD:2015]".to_string(), HashMap::from([("H", 20), ("C", 3), ("13C", 9), ("15N", 2), ("O", 2)])); // shTMT
    composition.insert("[UNIMOD:2016]".to_string(), HashMap::from([("H", 25), ("C", 8), ("13C", 7), ("N", 1), ("15N", 2), ("O", 3)])); // TMTpro
    composition.insert("[UNIMOD:2017]".to_string(), HashMap::from([("H", 25), ("C", 15), ("N", 3), ("O", 3)])); // TMTpro_zero
    composition.insert("[UNIMOD:2022]".to_string(), HashMap::from([("H", 12), ("C", 8), ("O", 7)])); // Ser/Thr-KDO
    composition.insert("[UNIMOD:2025]".to_string(), HashMap::from([("H", 28), ("C", 20), ("O", 4)])); // Andro-H2O
    composition.insert("[UNIMOD:2027]".to_string(), HashMap::from([("H", 7), ("C", 6), ("N", 3), ("O", 3)])); // His+O(2)
    composition.insert("[UNIMOD:2028]".to_string(), HashMap::from([("H", 176), ("C", 109), ("N", 8), ("O", 79)])); // Hex(6)HexNAc(5)NeuAc(3)
    composition.insert("[UNIMOD:2029]".to_string(), HashMap::from([("H", 148), ("C", 90), ("N", 6), ("O", 65)])); // Hex(7)HexNAc(6)
    composition.insert("[UNIMOD:2033]".to_string(), HashMap::from([("H", 9), ("C", 5), ("N", 1), ("O", 3), ("S", 1)])); // Met+O(2)
    composition.insert("[UNIMOD:2034]".to_string(), HashMap::from([("H", 3), ("C", 2), ("N", 1), ("O", 3)])); // Gly+O(2)
    composition.insert("[UNIMOD:2035]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 3)])); // Pro+O(2)
    composition.insert("[UNIMOD:2036]".to_string(), HashMap::from([("H", 12), ("C", 6), ("N", 2), ("O", 3)])); // Lys+O(2)
    composition.insert("[UNIMOD:2037]".to_string(), HashMap::from([("H", 7), ("C", 5), ("N", 1), ("O", 5)])); // Glu+O(2)
    composition.insert("[UNIMOD:2039]".to_string(), HashMap::from([("H", 24), ("C", 22), ("O", 8)])); // LTX+Lophotoxin
    composition.insert("[UNIMOD:2040]".to_string(), HashMap::from([("H", 108), ("C", 81), ("N", 7), ("O", 19)])); // MBS+peptide
    composition.insert("[UNIMOD:2041]".to_string(), HashMap::from([("H", 7), ("C", 7), ("O", 4), ("P", 1)])); // 3-hydroxybenzyl-phosphate
    composition.insert("[UNIMOD:2042]".to_string(), HashMap::from([("H", 5), ("C", 6), ("O", 3), ("P", 1)])); // phenyl-phosphate
    composition.insert("[UNIMOD:2044]".to_string(), HashMap::from([("H", 12), ("C", 9), ("N", 2), ("O", 6)])); // RBS-ID_Uridine
    composition.insert("[UNIMOD:2050]".to_string(), HashMap::from([("H", 25), ("13C", 15), ("15N", 3), ("O", 3)])); // shTMTpro
    composition.insert("[UNIMOD:2052]".to_string(), HashMap::from([("H", 70), ("C", 42), ("N", 8), ("O", 11), ("S", 1), ("Si", 1)])); // Biotin:Aha-DADPS
    composition.insert("[UNIMOD:2053]".to_string(), HashMap::from([("H", 38), ("C", 29), ("N", 8), ("O", 10), ("S", 1)])); // Biotin:Aha-PC
    composition.insert("[UNIMOD:2054]".to_string(), HashMap::from([("H", 10), ("C", 9), ("N", 2), ("O", 5)])); // pRBS-ID_4-thiouridine
    composition.insert("[UNIMOD:2055]".to_string(), HashMap::from([("H", 11), ("C", 10), ("N", 5), ("O", 4)])); // pRBS-ID_6-thioguanosine
    composition.insert("[UNIMOD:2057]".to_string(), HashMap::from([("H", 16), ("C", 8), ("N", 1), ("O", 4), ("P", 1)])); // 6C-CysPAT
    composition.insert("[UNIMOD:2058]".to_string(), HashMap::from([("H", 3), ("C", 8), ("O", 5), ("P", 1)])); // Xlink:DSPP[210]
    composition.insert("[UNIMOD:2059]".to_string(), HashMap::from([("H", 5), ("C", 8), ("O", 6), ("P", 1)])); // Xlink:DSPP[228]
    composition.insert("[UNIMOD:2060]".to_string(), HashMap::from([("H", 14), ("C", 12), ("N", 1), ("O", 8), ("P", 1)])); // Xlink:DSPP[331]
    composition.insert("[UNIMOD:2061]".to_string(), HashMap::from([("H", 5), ("C", 8), ("N", 1), ("O", 5), ("P", 1)])); // Xlink:DSPP[226]
    composition.insert("[UNIMOD:2062]".to_string(), HashMap::from([("H", 24), ("C", 14), ("N", 4), ("O", 3)])); // DBIA
    composition.insert("[UNIMOD:2067]".to_string(), HashMap::from([("H", 44), ("C", 26), ("N", 8), ("O", 8)])); // Mono_Ngamma-propargyl-L-Gln_desthiobiotin
    composition.insert("[UNIMOD:2068]".to_string(), HashMap::from([("H", 51), ("C", 31), ("N", 9), ("O", 10)])); // Di_L-Glu_Ngamma-propargyl-L-Gln_desthiobiotin
    composition.insert("[UNIMOD:2069]".to_string(), HashMap::from([("H", 52), ("C", 31), ("N", 10), ("O", 9)])); // Di_L-Gln_Ngamma-propargyl-L-Gln_desthiobiotin
    composition.insert("[UNIMOD:2070]".to_string(), HashMap::from([("H", 8), ("C", 5), ("N", 2), ("O", 2)])); // L-Gln
    composition.insert("[UNIMOD:2072]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 3)])); // Glyceroyl
    composition.insert("[UNIMOD:2073]".to_string(), HashMap::from([("H", 14), ("C", 13), ("N", 5), ("O", 6), ("P", 1)])); // N6pAMP
    composition.insert("[UNIMOD:2074]".to_string(), HashMap::from([("H", 21), ("C", 21), ("N", 5), ("O", 3)])); // DABCYL-C2-maleimide
    composition.insert("[UNIMOD:2079]".to_string(), HashMap::from([("H", 1), ("C", 6), ("N", 3), ("O", 3)])); // NBF
    composition.insert("[UNIMOD:2080]".to_string(), HashMap::from([("H", 12), ("C", 9), ("O", 3)])); // DCP
    composition.insert("[UNIMOD:2081]".to_string(), HashMap::from([("C", 2)])); // Ethynyl
    composition.insert("[UNIMOD:2082]".to_string(), HashMap::from([("H", 29), ("C", 18), ("N", 7), ("O", 8)])); // QQTGG
    composition.insert("[UNIMOD:2083]".to_string(), HashMap::from([("H", 26), ("C", 18), ("N", 6), ("O", 8)])); // Pyro-QQTGG
    composition.insert("[UNIMOD:2084]".to_string(), HashMap::from([("H", 27), ("C", 17), ("N", 7), ("O", 8)])); // NQTGG
    composition.insert("[UNIMOD:2085]".to_string(), HashMap::from([("H", 60), ("C", 41), ("N", 12), ("O", 15)])); // DVFQQQTGG
    composition.insert("[UNIMOD:2086]".to_string(), HashMap::from([("H", 11), ("C", 6), ("N", 1), ("O", 1)])); // iST-NHS specific cysteine modification
    composition.insert("[UNIMOD:2088]".to_string(), HashMap::from([("C", -2), ("13C", 2), ("N", -1), ("15N", 1)])); // Label:13C(2)15N(1)
    composition.insert("[UNIMOD:2106]".to_string(), HashMap::from([("H", 37), ("C", 21), ("N", 5), ("O", 6)])); // DPIA
    composition.insert("[UNIMOD:2107]".to_string(), HashMap::from([("H", 4), ("C", 4), ("O", 2)])); // Acetoacetyl
    composition.insert("[UNIMOD:2108]".to_string(), HashMap::from([("H", 8), ("C", 5), ("O", 1)])); // Isovaleryl
    composition.insert("[UNIMOD:2109]".to_string(), HashMap::from([("H", 8), ("C", 5), ("O", 1)])); // 2-methylbutyryl
    composition.insert("[UNIMOD:2110]".to_string(), HashMap::from([("H", 6), ("C", 5), ("O", 1)])); // Tiglyl
    composition.insert("[UNIMOD:2111]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 3)])); // 3-methylglutaryl
    composition.insert("[UNIMOD:2112]".to_string(), HashMap::from([("H", 6), ("C", 6), ("O", 3)])); // 3-methylglutaconyl
    composition.insert("[UNIMOD:2113]".to_string(), HashMap::from([("H", 8), ("C", 6), ("O", 4)])); // 3-hydroxy-3-methylglutaryl
    composition.insert("[UNIMOD:2114]".to_string(), HashMap::from([("H", 4), ("C", 3), ("O", 2)])); // Lactylation
    composition.insert("[UNIMOD:2115]".to_string(), HashMap::from([("H", 2), ("C", 3), ("O", 2)])); // Pyruvoyl
    composition.insert("[UNIMOD:2116]".to_string(), HashMap::from([("C", 2), ("O", 2)])); // Glyoxylyl
    composition.insert("[UNIMOD:2117]".to_string(), HashMap::from([("H", 6), ("C", 5), ("O", 4)])); // Itaconatyl
    composition.insert("[UNIMOD:2118]".to_string(), HashMap::from([("H", 4), ("C", 5), ("O", 3)])); // Itaconyl
    composition.insert("[UNIMOD:2119]".to_string(), HashMap::from([("H", 12), ("C", 7), ("N", 2), ("O", 2)])); // ValGly
    composition.insert("[UNIMOD:2120]".to_string(), HashMap::from([("H", 8), ("C", 5), ("O", 1)])); // Pentanoyl
    composition.insert("[UNIMOD:2121]".to_string(), HashMap::from([("H", 10), ("C", 6), ("O", 1)])); // Hexanoyl
    composition.insert("[UNIMOD:2122]".to_string(), HashMap::from([("H", 20), ("C", 2), ("13C", 10), ("N", -1), ("15N", 3), ("O", 2)])); // Label:13C(6)15N(2)+TMT6plex
    composition.insert("[UNIMOD:2123]".to_string(), HashMap::from([("H", 25), ("C", 2), ("13C", 13), ("N", -1), ("15N", 4), ("O", 3)])); // Label:13C(6)15N(2)+TMTpro
    composition.insert("[UNIMOD:2126]".to_string(), HashMap::from([("H", 8), ("C", 10), ("N", 4), ("S", 1)])); // 2PCA-triazole-ethanethiol

    composition
}
