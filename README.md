# CTRU failure-rate test

Run the unified test:

```powershell
python CTRU.py
```

`CTRU_failure.py` contains the CTRU Ring-3, CTRU-Light v1/v2/v3, and
CTRU-Prime 0715 estimators. CTRU and Prime support geometric and
Satterthwaite block variance methods, with optional E8 volume-threshold
correction. The default Light entry is v1.
