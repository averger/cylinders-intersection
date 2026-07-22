# Vérification mécanique de la théorie

Deux scripts prouvent les affirmations de [`../THEORY.md`](../THEORY.md) :

| Script | Rôle | Prérequis |
|---|---|---|
| `verify_theory.py` | 25 contrôles symboliques (SymPy) et numériques (NumPy) de toutes les équations | `pip install sympy numpy` |
| `verify_engine.py` | Validation croisée du moteur Rust contre la référence NumPy (< 2·10⁻¹³ mm) | serveur lancé sur `:8787`, `pip install numpy matplotlib` |

```bash
python -m venv venv && ./venv/bin/pip install sympy numpy matplotlib
./venv/bin/python verify_theory.py

cargo run --release -- --port 8787 &   # depuis la racine du dépôt
./venv/bin/python verify_engine.py
```

Chaque script sort avec un code non nul au moindre écart.
