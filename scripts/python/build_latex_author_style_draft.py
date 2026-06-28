"""Build a near-final LaTeX draft from the author-style manuscript.

The output is intended for expert feedback before journal submission. Citation
keys and BibTeX entries are generated from the local Mark-style bibliography,
but final public archive metadata should still be checked before submission.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_DIR = ROOT / "results" / "manuscript_package"
OUT_DIR = PACKAGE_DIR / "latex"
SOURCE_MD = PACKAGE_DIR / "manuscript_v1_author_style.md"
GW_PROJECTION_DIR = ROOT / "results" / "gw250114_constraints_comparison"


CITATION_REPLACEMENTS = [
    (
        "Black-hole ringdown spectroscopy is becoming an important test",
        r"Black-hole ringdown spectroscopy is becoming an important test~\cite{Berti2009QNMReview,Isi2019NoHair,Ghosh2021QNMConstraints}",
    ),
    (
        "GW250114 is a useful event for such a study because public spectroscopy\nproducts are available for more than one ringdown analysis",
        r"GW250114 is a useful event for such a study because public spectroscopy products are available for more than one ringdown analysis~\cite{LVK_GW250114_Spectroscopy,GW250114_ZenodoRelease}",
    ),
    (
        "higher-derivative Kerr QNM fingerprints",
        r"higher-derivative Kerr QNM fingerprints~\cite{Cano2024BeyondKerrQNM,Cano2023RotatingQNM,BeyondKerrQNMRepo}",
    ),
    (
        "Schwarzschild, braneworld tidal charge, Bardeen, and Hayward backgrounds",
        r"Schwarzschild, braneworld tidal charge, Bardeen, and Hayward backgrounds~\cite{Berti2009QNMReview,Toshmatov2016TidalCharge,Ulhoa2013Bardeen,Bolokhov2025Hayward}",
    ),
    (
        "For each mode, the Kerr reference is represented",
        r"For each mode, the Kerr reference is represented~\cite{Berti2006Spectroscopy,Stein2019QNM}",
    ),
    (
        "Many phenomenological black-hole\nmetrics are used in QPO, shadow, or geodesic studies",
        r"Many phenomenological black-hole metrics are used in QPO, shadow, or geodesic studies~\cite{Johannsen2015Metric,KRZ2016Parametrization,Pedrotti2024Eikonal}",
    ),
    (
        "Dynamical Chern-Simons or other theory-backed rotating calculations",
        r"Dynamical Chern-Simons or other theory-backed rotating calculations~\cite{Chung2025DCSMetrics,Li2025DCSSlowRotation,Pierini2022EdGB}",
    ),
]


BIBTEX = r"""
@article{Berti2009QNMReview,
  author = {Berti, Emanuele and Cardoso, Vitor and Starinets, Andrei O.},
  title = {Quasinormal modes of black holes and black branes},
  journal = {Classical and Quantum Gravity},
  volume = {26},
  pages = {163001},
  year = {2009},
  eprint = {0905.2975},
  archivePrefix = {arXiv},
  primaryClass = {gr-qc},
  note = {Verify final bibliographic fields before submission}
}

@article{Berti2006Spectroscopy,
  author = {Berti, Emanuele and Cardoso, Vitor and Will, Clifford M.},
  title = {On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA},
  journal = {Physical Review D},
  volume = {73},
  pages = {064030},
  year = {2006},
  eprint = {gr-qc/0512160},
  archivePrefix = {arXiv},
  note = {Used for Kerr QNM fitting-formula context}
}

@article{Stein2019QNM,
  author = {Stein, Leo C.},
  title = {{qnm}: A Python package for calculating Kerr quasinormal modes, separation constants, and spherical-spheroidal mixing coefficients},
  journal = {Journal of Open Source Software},
  volume = {4},
  pages = {1683},
  year = {2019},
  doi = {10.21105/joss.01683},
  note = {Verify version used locally}
}

@article{Isi2019NoHair,
  author = {Isi, Maximiliano and Giesler, Matthew and Farr, Will M. and Scheel, Mark A. and Teukolsky, Saul A.},
  title = {Testing the no-hair theorem with GW150914},
  journal = {Physical Review Letters},
  volume = {123},
  pages = {111102},
  year = {2019},
  eprint = {1905.00869},
  archivePrefix = {arXiv}
}

@article{Ghosh2021QNMConstraints,
  author = {Ghosh, Abhirup and others},
  title = {Constraints on quasinormal-mode frequencies with LIGO-Virgo binary-black-hole observations},
  year = {2021},
  eprint = {2104.01906},
  archivePrefix = {arXiv},
  note = {Verify journal and bibliographic details}
}

@article{LVK_GW250114_Spectroscopy,
  author = {{LVK Collaboration}},
  title = {Black Hole Spectroscopy and Tests of General Relativity with GW250114},
  year = {2026},
  note = {TODO: replace with final author list, arXiv identifier, journal, and DOI}
}

@misc{GW250114_ZenodoRelease,
  author = {{LVK Collaboration}},
  title = {GW250114 public posterior and spectroscopy data release},
  year = {2026},
  doi = {10.5281/zenodo.16877101},
  note = {Zenodo record 16877102; verify final citation text}
}

@article{Cano2024BeyondKerrQNM,
  author = {Cano, Pablo A. and others},
  title = {Complete higher-derivative corrections to Kerr quasinormal modes},
  year = {2024},
  eprint = {2409.04517},
  archivePrefix = {arXiv},
  note = {Primary QNM-fingerprint source; verify final journal fields}
}

@article{Cano2023RotatingQNM,
  author = {Cano, Pablo A. and others},
  title = {Rotating black-hole quasinormal modes in higher-derivative gravity},
  year = {2023},
  eprint = {2307.07431},
  archivePrefix = {arXiv},
  note = {Verify final journal fields}
}

@misc{BeyondKerrQNMRepo,
  author = {Cano, Pablo A. and collaborators},
  title = {{BeyondKerrQNM} public repository},
  year = {2024},
  howpublished = {\url{https://github.com/pacmn91/BeyondKerrQNM}},
  note = {Local import records commit 0afe6281bec6a6224bfd55fe4600d5966c6a7135}
}

@article{Toshmatov2016TidalCharge,
  author = {Toshmatov, Bobir and Stuchlik, Zdenek and Ahmedov, Bobomurat},
  title = {Quasinormal frequencies of black hole in the braneworld},
  year = {2016},
  eprint = {1605.02058},
  archivePrefix = {arXiv},
  note = {Verify final journal fields}
}

@article{Ulhoa2013Bardeen,
  author = {Ulhoa, S. C.},
  title = {On gravitational perturbations of regular black holes},
  year = {2013},
  eprint = {1303.3143},
  archivePrefix = {arXiv},
  note = {Bardeen axial gravitational WKB reference; verify final journal fields}
}

@article{MorenoSarbach2002NED,
  author = {Moreno, Claudia and Sarbach, Olivier},
  title = {Stability properties of black holes in self-gravitating nonlinear electrodynamics},
  year = {2002},
  eprint = {gr-qc/0208090},
  archivePrefix = {arXiv}
}

@article{Zhao2023Bardeen,
  author = {Zhao, Zheng and others},
  title = {Quasinormal modes of Bardeen black holes},
  year = {2023},
  eprint = {2306.02332},
  archivePrefix = {arXiv},
  note = {Verify title and journal fields}
}

@article{Bolokhov2025Hayward,
  author = {Bolokhov, S. V. and Skvortsova, M.},
  title = {Quasinormal modes of Hayward black holes},
  year = {2025},
  eprint = {2508.19989},
  archivePrefix = {arXiv},
  note = {Hayward WKB-Pade/time-domain reference; verify final journal fields}
}

@article{Pedraza2021Hayward,
  author = {Pedraza, O. and others},
  title = {Quasinormal modes of Hayward black holes with quintessence},
  year = {2021},
  eprint = {2111.06488},
  archivePrefix = {arXiv},
  note = {Verify final bibliographic fields}
}

@article{Johannsen2015Metric,
  author = {Johannsen, Tim},
  title = {Regular black hole metric with three constants of motion},
  year = {2015},
  eprint = {1501.02809},
  archivePrefix = {arXiv},
  note = {Metric/geodesic framework; not a gravitational-QNM source}
}

@article{KRZ2016Parametrization,
  author = {Konoplya, Roman and Rezzolla, Luciano and Zhidenko, Alexander},
  title = {General parametrization of axisymmetric black holes in metric theories of gravity},
  year = {2016},
  eprint = {1602.02378},
  archivePrefix = {arXiv}
}

@article{Pedrotti2024Eikonal,
  author = {Pedrotti, D. and Vagnozzi, S.},
  title = {Rotating regular black holes and the QNM-shadow correspondence},
  year = {2024},
  eprint = {2404.07589},
  archivePrefix = {arXiv},
  note = {Eikonal/geodesic reference; not used as observed 220/221 gravitational spectrum}
}

@article{Chung2025DCSMetrics,
  author = {Chung, A. K.-W. and Lam, H. and Yunes, N.},
  title = {Dynamical Chern-Simons quasinormal modes with METRICS},
  year = {2025},
  eprint = {2503.11759},
  archivePrefix = {arXiv},
  note = {Future theory-backed rotating extension}
}

@article{Li2025DCSSlowRotation,
  author = {Li, D. and Wagle, P. and Chen, H. and Yunes, N.},
  title = {Slow-rotation dynamical Chern-Simons quasinormal modes},
  year = {2025},
  eprint = {2503.15606},
  archivePrefix = {arXiv},
  note = {Future theory-backed overtone extension}
}

@article{Pierini2022EdGB,
  author = {Pierini, L. and Gualtieri, L.},
  title = {Quasinormal modes of rotating black holes in Einstein-dilaton-Gauss-Bonnet gravity},
  year = {2022},
  eprint = {2207.11267},
  archivePrefix = {arXiv},
  note = {Future theory-backed rotating extension}
}
""".strip()


USER_BIB_PATH = Path(
    r"C:\Users\vrb0015\Dropbox\Jarek\projekty\Magnetosphere of NS\Paper extrem dipoles\BH_mimickers_in_ph-sphere\references.bib"
)


CITATION_KEY_MAP = {
    "Berti2009QNMReview": "Ber-Car-Sta:2009:CQG:",
    "Berti2006Spectroscopy": "Ber-Car-Wil:2006:PRD:",
    "Stein2019QNM": "Ste:2019:JOSS:",
    "Isi2019NoHair": "Isi-Gie-Far:2019:PRL:",
    "Ghosh2021QNMConstraints": "Gho-Bri-Buo:2021:PRD:",
    "LVK_GW250114_Spectroscopy": "LIG-Vir-KAG:2025:ARXIV:",
    "GW250114_ZenodoRelease": "LIG-Vir-KAG:2025:ZENODO:",
    "Cano2024BeyondKerrQNM": "Can-Cap-Fra:2024:ARXIV:",
    "Cano2023RotatingQNM": "Can-Fra-Her:2023:ARXIV:",
    "BeyondKerrQNMRepo": "Can-etal:2024:GITHUB:",
    "Toshmatov2016TidalCharge": "Tos-Stu-Sch:2016:PRD:",
    "Ulhoa2013Bardeen": "Ulh:2014:BJP:",
    "MorenoSarbach2002NED": "Mor-Sar:2003:PRD:",
    "Zhao2023Bardeen": "Zha-etal:2023:ARXIV:",
    "Bolokhov2025Hayward": "Bol-Skv:2025:ARXIV:",
    "Pedraza2021Hayward": "Ped-Lop-Arc:2022:MPLA:",
    "Johannsen2015Metric": "Joh:2013:PRD:",
    "KRZ2016Parametrization": "Kon-Rez-Zhi:2016:PRD:",
    "Pedrotti2024Eikonal": "Ped-Vag:2024:PRD:",
    "Chung2025DCSMetrics": "Chu-Lam-Yun:2025:PRD:",
    "Li2025DCSSlowRotation": "Li-Wag-Che:2025:ARXIV:",
    "Pierini2022EdGB": "Pie-Gua:2022:PRD:",
}


MARK_STYLE_BIBTEX = r"""
@ARTICLE{Ber-Car-Sta:2009:CQG:,
       author = {{Berti}, Emanuele and {Cardoso}, Vitor and {Starinets}, Andrei O.},
        title = "{Quasinormal modes of black holes and black branes}",
      journal = {Classical and Quantum Gravity},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Phenomenology, High Energy Physics - Theory},
         year = 2009,
        month = aug,
       volume = {26},
       number = {16},
          eid = {163001},
        pages = {163001},
          doi = {10.1088/0264-9381/26/16/163001},
archivePrefix = {arXiv},
       eprint = {0905.2975},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2009CQGra..26p3001B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Ber-Car-Wil:2006:PRD:,
       author = {{Berti}, Emanuele and {Cardoso}, Vitor and {Will}, Clifford M.},
        title = "{On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA}",
      journal = {\prd},
     keywords = {04.70.-s, 04.30.Db, 04.80.Nn, Physics of black holes, Wave generation and sources, Gravitational wave detectors and experiments, General Relativity and Quantum Cosmology},
         year = 2006,
        month = mar,
       volume = {73},
       number = {6},
          eid = {064030},
        pages = {064030},
          doi = {10.1103/PhysRevD.73.064030},
archivePrefix = {arXiv},
       eprint = {gr-qc/0512160},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2006PhRvD..73f4030B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Ste:2019:JOSS:,
       author = {{Stein}, Leo C.},
        title = "{qnm: A Python package for calculating Kerr quasinormal modes, separation constants, and spherical-spheroidal mixing coefficients}",
      journal = {Journal of Open Source Software},
     keywords = {General Relativity and Quantum Cosmology, Instrumentation and Methods for Astrophysics},
         year = 2019,
        month = oct,
       volume = {4},
       number = {42},
          eid = {1683},
        pages = {1683},
          doi = {10.21105/joss.01683},
archivePrefix = {arXiv},
       eprint = {1908.10377},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019JOSS....4.1683S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Isi-Gie-Far:2019:PRL:,
       author = {{Isi}, Maximiliano and {Giesler}, Matthew and {Farr}, Will M. and {Scheel}, Mark A. and {Teukolsky}, Saul A.},
        title = "{Testing the no-hair theorem with GW150914}",
      journal = {\prl},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2019,
        month = sep,
       volume = {123},
       number = {11},
          eid = {111102},
        pages = {111102},
          doi = {10.1103/PhysRevLett.123.111102},
archivePrefix = {arXiv},
       eprint = {1905.00869},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019PhRvL.123k1102I},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Gho-Bri-Buo:2021:PRD:,
       author = {{Ghosh}, Abhirup and {Brito}, Richard and {Buonanno}, Alessandra},
        title = "{Constraints on quasi-normal-mode frequencies with LIGO-Virgo binary-black-hole observations}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2021,
        month = jun,
       volume = {103},
       number = {12},
          eid = {124041},
        pages = {124041},
          doi = {10.1103/PhysRevD.103.124041},
archivePrefix = {arXiv},
       eprint = {2104.01906},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021PhRvD.103l4041G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{LIG-Vir-KAG:2025:ARXIV:,
       author = {{LIGO Scientific Collaboration} and {Virgo Collaboration} and {KAGRA Collaboration}},
        title = "{Black Hole Spectroscopy and Tests of General Relativity with GW250114}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2025,
        month = sep,
          eid = {arXiv:2509.08099},
        pages = {arXiv:2509.08099},
archivePrefix = {arXiv},
       eprint = {2509.08099},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250908099L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@MISC{LIG-Vir-KAG:2025:ZENODO:,
       author = {{LIGO Scientific Collaboration} and {Virgo Collaboration} and {KAGRA Collaboration}},
        title = "{GW250114 discovery paper figure scripts}",
         year = 2025,
        month = sep,
          doi = {10.5281/zenodo.16877102},
 howpublished = {Zenodo},
         note = {Version v1; contains code and data required to reproduce the figures in the GW250114 discovery paper and supplement}
}

@ARTICLE{Can-Cap-Fra:2024:ARXIV:,
       author = {{Cano}, Pablo A. and {Capuano}, Lodovico and {Franchini}, Nicola and {Maenaut}, Simon and {V{\"o}lkel}, Sebastian H.},
        title = "{Higher-derivative corrections to the Kerr quasinormal mode spectrum}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Theory},
         year = 2024,
        month = sep,
          eid = {arXiv:2409.04517},
        pages = {arXiv:2409.04517},
archivePrefix = {arXiv},
       eprint = {2409.04517},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv240904517C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Can-Fra-Her:2023:ARXIV:,
       author = {{Cano}, Pablo A. and {Fransen}, Kwinten and {Hertog}, Thomas and {Maenaut}, Simon},
        title = "{Quasinormal modes of rotating black holes in higher-derivative gravity}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Theory},
         year = 2023,
        month = jul,
          eid = {arXiv:2307.07431},
        pages = {arXiv:2307.07431},
archivePrefix = {arXiv},
       eprint = {2307.07431},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230707431C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@MISC{Can-etal:2024:GITHUB:,
       author = {{Cano}, Pablo A. and collaborators},
        title = "{BeyondKerrQNM public repository}",
         year = 2024,
 howpublished = {\url{https://github.com/pacmn91/BeyondKerrQNM}},
         note = {Local import records commit 0afe6281bec6a6224bfd55fe4600d5966c6a7135}
}

@ARTICLE{Joh:2013:PRD:,
       author = {{Johannsen}, Tim},
        title = "{Regular Black Hole Metric with Three Constants of Motion}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2013,
        month = aug,
       volume = {88},
       number = {4},
          eid = {044002},
        pages = {044002},
          doi = {10.1103/PhysRevD.88.044002},
archivePrefix = {arXiv},
       eprint = {1501.02809},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2013PhRvD..88d4002J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Kon-Rez-Zhi:2016:PRD:,
       author = {{Konoplya}, Roman and {Rezzolla}, Luciano and {Zhidenko}, Alexander},
        title = "{General parametrization of axisymmetric black holes in metric theories of gravity}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Theory},
         year = 2016,
        month = mar,
       volume = {93},
       number = {6},
          eid = {064015},
        pages = {064015},
          doi = {10.1103/PhysRevD.93.064015},
archivePrefix = {arXiv},
       eprint = {1602.02378},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016PhRvD..93f4015K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Ped-Vag:2024:PRD:,
       author = {{Pedrotti}, Davide and {Vagnozzi}, Sunny},
        title = "{Quasinormal modes-shadow correspondence for rotating regular black holes}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Phenomenology, High Energy Physics - Theory},
         year = 2024,
        month = oct,
       volume = {110},
       number = {8},
          eid = {084075},
        pages = {084075},
          doi = {10.1103/PhysRevD.110.084075},
archivePrefix = {arXiv},
       eprint = {2404.07589},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024PhRvD.110h4075P},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Tos-Stu-Sch:2016:PRD:,
       author = {{Toshmatov}, Bobir and {Stuchl{\'i}k}, Zden{\v e}k and {Schee}, Jan and {Ahmedov}, Bobomurat},
        title = "{Quasinormal frequencies of black hole in the braneworld}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology},
         year = 2016,
        month = jun,
       volume = {93},
       number = {12},
          eid = {124017},
        pages = {124017},
          doi = {10.1103/PhysRevD.93.124017},
archivePrefix = {arXiv},
       eprint = {1605.02058},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016PhRvD..93l4017T},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Ulh:2014:BJP:,
       author = {{Ulhoa}, S. C.},
        title = "{On Quasinormal Modes for Gravitational Perturbations of Bardeen Black Hole}",
      journal = {Brazilian Journal of Physics},
     keywords = {General Relativity and Quantum Cosmology},
         year = 2014,
        month = aug,
       volume = {44},
       number = {4},
        pages = {380-384},
          doi = {10.1007/s13538-014-0209-7},
archivePrefix = {arXiv},
       eprint = {1303.3143},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014BrJPh..44..380U},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Bol-Skv:2025:ARXIV:,
       author = {{Bolokhov}, S. V. and {Skvortsova}, Milena},
        title = "{Gravitational quasinormal modes of the Hayward spacetime}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology},
         year = 2025,
        month = aug,
          eid = {arXiv:2508.19989},
        pages = {arXiv:2508.19989},
archivePrefix = {arXiv},
       eprint = {2508.19989},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250819989B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Chu-Lam-Yun:2025:PRD:,
       author = {{Chung}, Adrian Ka-Wai and {Lam}, Kelvin Ka-Ho and {Yunes}, Nicolas},
        title = "{Quasinormal mode frequencies and gravitational perturbations of spinning black holes in modified gravity through METRICS: The dynamical Chern-Simons gravity case}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, Instrumentation and Methods for Astrophysics},
         year = 2025,
        month = jun,
       volume = {111},
       number = {12},
          eid = {124052},
        pages = {124052},
          doi = {10.1103/g83n-rrlj},
archivePrefix = {arXiv},
       eprint = {2503.11759},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025PhRvD.111l4052C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Li-Wag-Che:2025:ARXIV:,
       author = {{Li}, Dongjun and {Wagle}, Pratik and {Chen}, Yanbei and {Yunes}, Nicol{\'a}s},
        title = "{Perturbations of spinning black holes in dynamical Chern-Simons gravity: Slow rotation quasinormal modes}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology},
         year = 2025,
        month = mar,
          eid = {arXiv:2503.15606},
        pages = {arXiv:2503.15606},
archivePrefix = {arXiv},
       eprint = {2503.15606},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250315606L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Pie-Gua:2022:PRD:,
       author = {{Pierini}, Lorenzo and {Gualtieri}, Leonardo},
        title = "{Quasi-normal modes of rotating black holes in Einstein-dilaton Gauss-Bonnet gravity: the second order in rotation}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, High Energy Physics - Theory},
         year = 2022,
        month = nov,
       volume = {106},
       number = {10},
          eid = {104009},
        pages = {104009},
          doi = {10.1103/PhysRevD.106.104009},
archivePrefix = {arXiv},
       eprint = {2207.11267},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022PhRvD.106j4009P},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Mor-Sar:2003:PRD:,
       author = {{Moreno}, Claudia and {Sarbach}, Olivier},
        title = "{Stability properties of black holes in self-gravitating nonlinear electrodynamics}",
      journal = {\prd},
     keywords = {General Relativity and Quantum Cosmology},
         year = 2003,
        month = jan,
       volume = {67},
       number = {2},
          eid = {024028},
        pages = {024028},
          doi = {10.1103/PhysRevD.67.024028},
archivePrefix = {arXiv},
       eprint = {gr-qc/0208090},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2003PhRvD..67b4028M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{Ped-Lop-Arc:2022:MPLA:,
       author = {{Pedraza}, Omar and {L{\'o}pez}, L. A. and {Arceo}, R. and {Cabrera-Munguia}, I.},
        title = "{Quasinormal modes of the Hayward black hole surrounded by quintessence: scalar, electromagnetic and gravitational perturbations}",
      journal = {Modern Physics Letters A},
     keywords = {General Relativity and Quantum Cosmology, High Energy Physics - Phenomenology},
         year = 2022,
          doi = {10.1142/S0217732322500572},
archivePrefix = {arXiv},
       eprint = {2111.06488},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MPLA...3750057P},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
""".strip()


CITATION_MAP = """# Citation Map For The LaTeX Draft

Generated for `manuscript_v1_author_style.tex`.

## Core ringdown and Kerr references

- `Berti2009QNMReview`: broad QNM/ringdown review and notation background.
- `Berti2006Spectroscopy`: Kerr spectroscopy fitting-formula context and Berti table usage.
- `Stein2019QNM`: numerical Kerr QNM cross-check through the `qnm` package.
- `Isi2019NoHair`, `Ghosh2021QNMConstraints`: observational ringdown spectroscopy and QNM-constraint context.

## Public GW250114 inputs

- `LVK_GW250114_Spectroscopy`: main GW250114 spectroscopy/results paper; bibliographic fields must be replaced by final public metadata.
- `GW250114_ZenodoRelease`: public posterior/spectroscopy products used by the projection branch.

## Higher-derivative Kerr branch

- `Cano2024BeyondKerrQNM`, `Cano2023RotatingQNM`: QNM fingerprints used for the public-data projection.
- `BeyondKerrQNMRepo`: public fit-data/code source; commit and data-import record should be checked against the local provenance files.

## Static supplied-potential branch

- `Toshmatov2016TidalCharge`: braneworld tidal-charge QNM benchmark.
- `Ulhoa2013Bardeen`, `MorenoSarbach2002NED`, `Zhao2023Bardeen`: Bardeen/NED perturbation and QNM context.
- `Bolokhov2025Hayward`, `Pedraza2021Hayward`: Hayward QNM/overtone and comparison context.

## Phenomenological metric motivation

- `Johannsen2015Metric`, `KRZ2016Parametrization`, `Pedrotti2024Eikonal`: examples of metric/geodesic or eikonal motivation; these are not used as gravitational 220/221 spectra.

## Outlook candidates

- `Chung2025DCSMetrics`, `Li2025DCSSlowRotation`, `Pierini2022EdGB`: theory-backed rotating QNM directions for future extensions.

## Verification still needed

- Replace all `TODO` BibTeX notes with final journal/DOI/arXiv metadata.
- Check exact title and author spelling for every static-metric reference.
- Add citations for every table row once the final table captions are fixed.
- Decide whether PRD/PRL-style `apsrev4-2` is the final bibliography style.
"""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def latex_escape(text: object) -> str:
    s = str(text)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    return s


def fmt_float(value: object, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return latex_escape(value)


def operator_label(operator: str, polarization: str | None = None) -> str:
    labels = {
        "lambda_ev": r"\lambda_{\rm ev}",
        "lambda_odd": r"\lambda_{\rm odd}",
        "epsilon1": r"\epsilon_1",
        "epsilon2": r"\epsilon_2",
        "epsilon3": r"\epsilon_3",
    }
    label = labels.get(operator, latex_escape(operator))
    if polarization:
        signs = {"plus": "+", "minus": "-"}
        sign = signs.get(polarization, latex_escape(polarization))
        return rf"${label}^{{{sign}}}$"
    return rf"${label}$"


def projection_label(projection: str) -> str:
    labels = {
        "RINGDOWN": "RINGDOWN",
        "PYRING_DELTA": "pyRing deviation",
        "RINGDOWN_vs_PYRING": "RINGDOWN vs. pyRing",
    }
    return labels.get(projection, latex_escape(projection))


def table1_latex() -> str:
    rows = read_csv(PACKAGE_DIR / "main_table1_public_projection_summary.csv")
    observable_labels = {
        "{log f_220, log f_221, df_221}": r"$\log f_{220}$, $\log f_{221}$, $\delta f_{221}$",
        "{log(1+domega_221), log(1+dtau_221)}": r"$\delta\omega_{221}$, $\delta\tau_{221}$",
        "operator-by-operator comparison": "operator comparison",
    }
    row_labels = {
        "lambda_odd plus": r"$\lambda_{\rm odd}^{+}$",
        "epsilon1 plus": r"$\epsilon_1^{+}$",
    }
    body = []
    for row in rows:
        body.append(
            " & ".join(
                [
                    projection_label(row["projection"]),
                    observable_labels.get(row["observable_set"], latex_escape(row["observable_set"])),
                    latex_escape(row["tested_rows"]),
                    latex_escape(row["zero_outside_90pct"]),
                    latex_escape(row["max_sigma_from_zero"]),
                    row_labels.get(row["max_sigma_row"], latex_escape(row["max_sigma_row"])),
                ]
            )
            + r" \\"
        )
    return "\n".join(
        [
            r"\begin{table*}",
            r"\caption{Public GW250114 higher-derivative QNM projections. The RINGDOWN and pyRing branches are reported separately and are compared but not statistically combined.}",
            r"\label{tab:public_projection_summary}",
            r"\begin{ruledtabular}",
            r"\begin{tabular}{lccccc}",
            r"Projection & Observables & Rows & Zero outside 90\% & Max statistic & Row \\",
            *body,
            r"\end{tabular}",
            r"\end{ruledtabular}",
            r"\end{table*}",
        ]
    )


def projection_detail_latex() -> str:
    rows = read_csv(GW_PROJECTION_DIR / "projection_constraints_long.csv")
    body = []
    for row in rows:
        branch = projection_label(row["projection"])
        op = operator_label(row["operator"], row["polarization"])
        alpha_best = fmt_float(row["alpha_best"], 3)
        alpha_sigma = fmt_float(row["alpha_sigma"], 3)
        interval = rf"$[{fmt_float(row['alpha_q05'], 3)},\,{fmt_float(row['alpha_q95'], 3)}]$"
        sigma_from_zero = fmt_float(row["sigma_from_alpha0"], 3)
        zero_inside = "yes" if row["zero_inside_90pct"] == "1" else "no"
        body.append(
            " & ".join([branch, op, alpha_best, alpha_sigma, interval, sigma_from_zero, zero_inside])
            + r" \\"
        )
    return "\n".join(
        [
            r"\begin{table*}",
            r"\caption{Per-operator public GW250114 projection values. The intervals are one-at-a-time 90 percent public-product intervals for the local linearized coupling $\alpha$; they are not Bayes factors and not multi-parameter EFT constraints.}",
            r"\label{tab:public_projection_detail}",
            r"\begin{ruledtabular}",
            r"\begin{tabular}{llccccc}",
            r"Branch & Operator & $\hat{\alpha}$ & $\sigma_{\alpha}$ & 90\% interval & $|\hat{\alpha}|/\sigma_{\alpha}$ & $\alpha=0$ inside \\",
            *body,
            r"\end{tabular}",
            r"\end{ruledtabular}",
            r"\end{table*}",
        ]
    )


def projection_forest_figure_latex() -> str:
    return "\n".join(
        [
            r"\begin{figure*}[t]",
            r"\centering",
            r"\includegraphics[width=0.95\textwidth]{../../gw250114_constraints_comparison/projection_alpha_interval_forest.pdf}",
            r"\caption{Forest plot of the one-at-a-time public GW250114 projections. Points show $\hat{\alpha}$ and bars show the 90 percent public-product intervals. The vertical line marks $\alpha=0$, which lies inside every interval in both the RINGDOWN and pyRing branches.}",
            r"\label{fig:public_projection_forest}",
            r"\end{figure*}",
        ]
    )


def readiness_ladder_figure_latex() -> str:
    return "\n".join(
        [
            r"\begin{figure*}[t]",
            r"\centering",
            r"\IfFileExists{../../static_qnm_readiness_audit/static_metric_readiness_ladder.pdf}{%",
            r"\includegraphics[width=0.9\textwidth]{../../static_qnm_readiness_audit/static_metric_readiness_ladder.pdf}%",
            r"}{%",
            r"\fbox{\begin{minipage}{0.86\textwidth}",
            r"\centering Static readiness ladder available as",
            r"\path{../../static_qnm_readiness_audit/static_metric_readiness_ladder.svg}.",
            r"Convert the SVG to PDF before final journal submission.",
            r"\end{minipage}}%",
            r"}",
            r"\caption{Ringdown-readiness ladder for metric models. A line element or geodesic observable is not sufficient for gravitational ringdown use without perturbation physics and a validated spectrum.}",
            r"\label{fig:readiness_ladder}",
            r"\end{figure*}",
        ]
    )


def table2_latex() -> str:
    rows = read_csv(PACKAGE_DIR / "main_table2_static_readiness_summary.csv")
    family_labels = {
        "Bardeen_NED": "Bardeen NED",
        "Braneworld_tidal_charge": "Braneworld tidal charge",
        "Hayward": "Hayward",
        "Hayward_overtone": "Hayward overtone",
        "Schwarzschild": "Schwarzschild",
    }
    mode_labels = {
        "l=2,n=0": r"$\ell=2,n=0$",
        "l=2,n=1": r"$\ell=2,n=1$",
    }
    parameter_labels = {
        "alpha=0.7698": r"$\alpha=0.7698$",
        "q_tidal=2": r"$q_{\rm tidal}=2$",
        "gamma=1.18": r"$\gamma=1.18$",
        "not_applicable": "n/a",
        "": "n/a",
    }
    body = []
    for row in rows:
        physical = row["max_sampled_physical_delta_pct"] or "n/a"
        parameter = parameter_labels.get(row["largest_sampled_parameter"], latex_escape(row["largest_sampled_parameter"]))
        body.append(
            " & ".join(
                [
                    family_labels.get(row["family"], latex_escape(row["family"])),
                    mode_labels.get(row["validated_modes"], latex_escape(row["validated_modes"])),
                    latex_escape(row["validation_rows"]),
                    latex_escape(row["max_validation_delta_pct"]),
                    latex_escape(physical),
                    parameter,
                ]
            )
            + r" \\"
        )
    return "\n".join(
        [
            r"\begin{table*}",
            r"\caption{Static supplied-potential QNM readiness summary. Validation errors quantify reproduction of reference spectra. The maximum physical shift is the largest sampled absolute relative change in either $\mathrm{Re}(M\omega)$ or $|\mathrm{Im}(M\omega)|$ from the Schwarzschild or zero-parameter baseline; it is not an observational exclusion interval.}",
            r"\label{tab:static_readiness_summary}",
            r"\begin{ruledtabular}",
            r"\begin{tabular}{lccccc}",
            r"Family & Modes & Rows & Max validation [\%] & Max physical shift [\%] & Largest sampled parameter \\",
            *body,
            r"\end{tabular}",
            r"\end{ruledtabular}",
            r"\end{table*}",
        ]
    )


def apply_citations(md: str) -> str:
    return md


def apply_citation_key_map(text: str) -> str:
    out = text
    for old, new in CITATION_KEY_MAP.items():
        out = out.replace(old, new)
    return out


def project_references_bib() -> str:
    if USER_BIB_PATH.exists():
        base = USER_BIB_PATH.read_text(encoding="utf-8", errors="replace").rstrip()
        return "\n\n".join(
            [
                base,
                "% --- Added for GW250114/QNM-readiness project; keys follow the local Mark citation style. ---",
                MARK_STYLE_BIBTEX,
            ]
        )
    return MARK_STYLE_BIBTEX


def convert_inline(text: str) -> str:
    inline_map = {
        "alpha = 0": r"$\alpha=0$",
        "alpha": r"$\alpha$",
        "0.497 sigma": r"$0.497\sigma$",
        "1.222 sigma": r"$1.222\sigma$",
        "0.535 sigma": r"$0.535\sigma$",
        "1.290 sigma": r"$1.290\sigma$",
        "1.286": r"$1.286$",
        "25": r"$25$",
        "0.925%": r"$0.925\%$",
        "abs(median)/sd": r"$|\mathrm{median}|/\sigma$",
        "lambda_odd plus": r"$\lambda_{\rm odd}^{+}$",
        "epsilon1 plus": r"$\epsilon_1^{+}$",
        "alpha_hat": r"$\hat{\alpha}$",
        "sigma_alpha": r"$\sigma_{\alpha}$",
        "https://github.com/jvrba/gw250114-ringdown-readiness": r"\url{https://github.com/jvrba/gw250114-ringdown-readiness}",
        "y": r"$y$",
        "mu": r"$\mu$",
        "C": r"$C$",
        "theta = {delta log M_f, delta chi_f}": r"$\theta=\{\delta\log M_f,\delta\chi_f\}$",
        "s": r"$s$",
        "{-1, -1, 0}": r"$\{-1,-1,0\}$",
        "{log f_220, log f_221, df_221}": r"$\{\log f_{220},\log f_{221},\delta f_{221}\}$",
        "x": r"$x$",
        "qnm": r"\texttt{qnm}",
    }

    def replace_code(match: re.Match[str]) -> str:
        code = match.group(1)
        if code in inline_map:
            return inline_map[code]
        return r"\texttt{" + latex_escape(code) + "}"

    return re.sub(r"`([^`]+)`", replace_code, text)


def latex_block(content: str) -> list[str]:
    stripped = content.strip()
    math_blocks = {
        "y_R = {log f_220, log f_221, df_221}.": [
            r"\begin{equation}",
            r"y_R=\{\log f_{220},\log f_{221},\delta f_{221}\}.",
            r"\end{equation}",
        ],
        "y_P = {log(1 + domega_221), log(1 + dtau_221)}.": [
            r"\begin{equation}",
            r"y_P=\{\log(1+\delta\omega_{221}),\log(1+\delta\tau_{221})\}.",
            r"\end{equation}",
        ],
        "M omega_lmn(chi) = M omega_R,lmn(chi) - i M/tau_lmn(chi).": [
            r"\begin{equation}",
            r"M\omega_{\ell mn}(\chi)=M\omega_{R,\ell mn}(\chi)-i\,M/\tau_{\ell mn}(\chi).",
            r"\end{equation}",
        ],
        "f_lmn(M, chi) propto Re[M omega_lmn(chi)] / M,\ntau_lmn(M, chi) propto M / Abs(Im[M omega_lmn(chi)]).": [
            r"\begin{align}",
            r"f_{\ell mn}(M,\chi)&\propto \frac{\mathrm{Re}[M\omega_{\ell mn}(\chi)]}{M},\\",
            r"\tau_{\ell mn}(M,\chi)&\propto \frac{M}{|\mathrm{Im}[M\omega_{\ell mn}(\chi)]|}.",
            r"\end{align}",
        ],
        "d log f_lmn / d alpha,\nd log tau_lmn / d alpha.": [
            r"\begin{align}",
            r"\frac{\partial \log f_{\ell mn}}{\partial \alpha},\qquad",
            r"\frac{\partial \log \tau_{\ell mn}}{\partial \alpha}.",
            r"\end{align}",
        ],
        "chi2(p) = (mu - y_model(p))^T C^-1 (mu - y_model(p)).": [
            r"\begin{equation}",
            r"\chi^2(p)=\left[\mu-y_{\rm model}(p)\right]^T C^{-1}\left[\mu-y_{\rm model}(p)\right].",
            r"\end{equation}",
        ],
        "r = mu - y_Kerr = N theta + s alpha,": [
            r"\begin{equation}",
            r"r=\mu-y_{\rm Kerr}=N\theta+s\alpha,",
            r"\end{equation}",
        ],
        "y_P = s alpha,\ns = {d log f_221 / d alpha, d log tau_221 / d alpha}.": [
            r"\begin{align}",
            r"y_P&=s\alpha,\\",
            r"s&=\left\{\frac{\partial\log f_{221}}{\partial\alpha},",
            r"\frac{\partial\log\tau_{221}}{\partial\alpha}\right\}.",
            r"\end{align}",
        ],
        "alpha_hat +/- 1.6448536269514722 sigma_alpha.": [
            r"\begin{equation}",
            r"\hat{\alpha}\pm 1.6448536269514722\,\sigma_{\alpha}.",
            r"\end{equation}",
        ],
        "D = Abs(alpha_R - alpha_P) / Sqrt(sigma_R^2 + sigma_P^2).": [
            r"\begin{equation}",
            r"D=\frac{|\alpha_R-\alpha_P|}{\sqrt{\sigma_R^2+\sigma_P^2}}.",
            r"\end{equation}",
        ],
        "d_t^2 Psi - d_x^2 Psi + V(r(x)) Psi = 0,": [
            r"\begin{equation}",
            r"\partial_t^2\Psi-\partial_x^2\Psi+V(r(x))\Psi=0,",
            r"\end{equation}",
        ],
    }
    if stripped in math_blocks:
        return math_blocks[stripped]

    if stripped == "\n".join(
        [
            "metric/geodesic only",
            "test-field QNM",
            "supplied gravitational master potential",
            "validated static gravitational QNM",
            "theory-backed rotating gravitational QNM.",
        ]
    ):
        return [
            r"\begin{itemize}",
            r"\item metric/geodesic only;",
            r"\item test-field QNM;",
            r"\item supplied gravitational master potential;",
            r"\item validated static gravitational QNM;",
            r"\item theory-backed rotating gravitational QNM.",
            r"\end{itemize}",
        ]

    if stripped.startswith("Bardeen NED:"):
        return [
            r"\begin{itemize}",
            r"\item Bardeen NED: maximum sampled shift $109.8\%$ at $\alpha=0.7698$;",
            r"\item braneworld tidal charge: maximum sampled shift $21.214\%$ at $q_{\rm tidal}=2$;",
            r"\item Hayward fundamental: maximum sampled shift $16.596\%$ at $\gamma=1.18$;",
            r"\item Hayward first overtone: maximum sampled shift $17.407\%$ at $\gamma=1.18$.",
            r"\end{itemize}",
        ]

    return [
        r"\begin{verbatim}",
        *content.splitlines(),
        r"\end{verbatim}",
    ]


def md_to_latex(md: str) -> tuple[str, str]:
    lines = md.splitlines()
    title = lines[0].lstrip("# ").strip()
    body_lines = []
    in_code = False
    code_lines: list[str] = []
    skip_abstract_heading = False

    for raw in lines[1:]:
        line = raw.rstrip()
        if line.startswith("```"):
            if not in_code:
                in_code = True
                code_lines = []
            else:
                in_code = False
                content = "\n".join(code_lines)
                body_lines.extend(latex_block(content))
            continue

        if in_code:
            code_lines.append(line)
            continue

        if line == "## Abstract":
            skip_abstract_heading = True
            body_lines.append(r"\begin{abstract}")
            continue
        if line.startswith("## 1. Introduction") and skip_abstract_heading:
            body_lines.append(r"\end{abstract}")
            body_lines.append(r"\maketitle")
            body_lines.append(r"\section{Introduction}")
            continue
        if line.startswith("## "):
            heading = re.sub(r"^##\s*(\d+\.\s*)?", "", line)
            if heading in {"Data and Code Availability", "Data And Code Availability Draft", "Acknowledgement And Citation TODO"}:
                body_lines.append(r"\section*{" + latex_escape(heading.replace(" Draft", "")) + "}")
            else:
                body_lines.append(r"\section{" + latex_escape(heading) + "}")
            continue
        if line.startswith("### "):
            heading = re.sub(r"^###\s*(\d+\.\d+\s*)?", "", line)
            body_lines.append(r"\subsection{" + latex_escape(heading) + "}")
            continue
        if not line.strip():
            body_lines.append("")
            continue

        body_lines.append(convert_inline(line))

    return title, "\n".join(body_lines)


def build_tex() -> str:
    md = SOURCE_MD.read_text(encoding="utf-8")
    md = apply_citations(md)
    title, body = md_to_latex(md)
    body = body.replace(
        "rotating calculation.\n\n" + r"\subsection{Public GW250114 Inputs}",
        "rotating calculation.\n\n"
        + readiness_ladder_figure_latex()
        + "\n\n"
        + r"\subsection{Public GW250114 Inputs}",
    )
    body = body.replace(
        r"\subsection{Public GW250114 Projection}",
        projection_forest_figure_latex()
        + "\n\n"
        + table1_latex()
        + "\n\n"
        + projection_detail_latex()
        + "\n\n"
        + r"\subsection{Public GW250114 Projection}",
    )
    body = body.replace(
        r"\subsection{Static QNM Benchmarks}",
        table2_latex() + "\n\n" + r"\subsection{Static QNM Benchmarks}",
    )
    body = body.replace(r"\section{Results}", r"\section{Results and Discussion}")
    body = body.replace(r"\section{Discussion}", r"\subsection{Interpretation and Limitations}")
    tex = "\n".join(
        [
            r"\documentclass[aps,prd,10pt,twocolumn,nofootinbib,superscriptaddress]{revtex4-2}",
            r"\usepackage{amsmath,amssymb}",
            r"\usepackage{graphicx}",
            r"\usepackage{url}",
            r"\usepackage{hyperref}",
            "",
            r"\begin{document}",
            "",
            r"\title{" + latex_escape(title) + "}",
            r"\author{Author list to be completed}",
            r"\affiliation{Affiliation to be completed}",
            "",
            body,
            "",
            r"\clearpage",
            r"\bibliographystyle{apsrev4-2}",
            r"\bibliography{references_draft}",
            "",
            r"\end{document}",
            "",
        ]
    )
    return apply_citation_key_map(tex)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tex_path = OUT_DIR / "manuscript_v1_author_style.tex"
    bib_path = OUT_DIR / "references_draft.bib"
    notes_path = OUT_DIR / "latex_build_notes.md"
    citation_map_path = OUT_DIR / "citation_map.md"

    tex_path.write_text(build_tex(), encoding="utf-8")
    bib_path.write_text(project_references_bib() + "\n", encoding="utf-8")
    citation_map_path.write_text(apply_citation_key_map(CITATION_MAP) + "\n", encoding="utf-8")
    notes_path.write_text(
        "\n".join(
            [
                "# LaTeX Draft Notes",
                "",
                "Generated from `results/manuscript_package/manuscript_v1_author_style.md`.",
                "",
                "## Important",
                "",
                "- This is a near-final PRD-style expert-feedback draft.",
                "- The repository URL is still a placeholder and must be replaced by the final public archive before journal submission.",
                "- Bibliographic metadata in `references_draft.bib` should be checked one final time against ADS/journal records.",
                "- Author list, affiliations, acknowledgements, and exact repository/archive metadata are still the main remaining submission items.",
                "- Use `citation_map.md` for the final reference-verification pass.",
                "- The language follows the author-style profile; do not over-polish away the author's voice.",
                "",
                "## Suggested Next Commands",
                "",
                "```powershell",
                "cd results/manuscript_package/latex",
                "pdflatex manuscript_v1_author_style.tex",
                "bibtex manuscript_v1_author_style",
                "pdflatex manuscript_v1_author_style.tex",
                "pdflatex manuscript_v1_author_style.tex",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"LaTeX draft: {tex_path}")
    print(f"BibTeX draft: {bib_path}")
    print(f"Citation map: {citation_map_path}")
    print(f"Notes: {notes_path}")


if __name__ == "__main__":
    main()
