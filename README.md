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
