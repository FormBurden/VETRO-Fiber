# modules/nap/constants.py

# Tunables (NAP-focused)
VAULT_ON_LINE_FT: float = 15.0         # snap vaults to nearest conduit if within this distance
NODE_SNAP_FT: float = 20.0              # treat near-mid-segment tees as shared nodes if within this
NAP_MAX_SL: int = 4                     # max service locations per NAP

# While drops have their own distance limit, NAP logic uses the above exclusively.
