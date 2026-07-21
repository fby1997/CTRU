# CTRU failure-rate test

Run the unified test:

```powershell
python CTRU.py
```

`CTRU_failure.py` contains the CTRU, CTRU-Light, and
CTRU-Prime estimators. CTRU and Prime support geometric and
Satterthwaite block variance methods, with optional E8 volume-threshold
correction. The default Light entry is v1.

## Repository Updates & Notices

* **Branch Structure Update:** Moving forward, the `main` branch will exclusively host the test code for the latest version. All previous versions will be archived and maintained in their respective version-specific branches. 
* **CTRU Error Rate Calculation:** Please note that the current error rate calculation code for CTRU has been corrected and updated compared to our previously published papers. A detailed explanation and comprehensive analysis of these corrections will be provided in our upcoming paper.
**3. Updates to Security Analysis**
Furthermore, there may be variations in the security analysis results, most notably in the quantum security estimations. The security evaluations for the current version have been updated based on the findings from:
> Xavier Bonnetain, André Chailloux, André Schrottenloher, Yixin Shen. "[Finding many Collisions via Reusable Quantum Walks](https://doi.org/10.1007/978-3-031-30589-4_8)." *EUROCRYPT 2023*. [DOI: 10.1007/978-3-031-30589-4_8](https://doi.org/10.1007/978-3-031-30589-4_8) | ⟨[hal-04261002](https://hal.science/hal-04261002)⟩
