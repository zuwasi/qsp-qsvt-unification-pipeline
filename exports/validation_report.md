# QSP / QSVT Pipeline -- Validation Report

Paper: Grand Unification of Quantum Algorithms (Martyn, Rossi, Tan, Chuang, 2021).

| Check | Result |
|---|---|
| Chebyshev symbolic derivation (d=1..6) all pass | True |
| Zero-phase QSP vs T_d, max error |           -16
7.77156 10 |
| QSP unitarity |P|^2+(1-a^2)|Q|^2=1, max error |           -16
8.88178 10 |
| Jacobi-Anger cos truncation error |           -11
2.15186 10 |
| Jacobi-Anger sin truncation error |           -12
1.42161 10 |
| Phase generation, cos target max error |           -7
2.78753 10 |
| Phase generation, sin target max error |           -8
6.59942 10 |
| e^{-iHt} reconstruction max error |           -7
1.34261 10 |
| e^{-iHt} normalized overlap | 1. |