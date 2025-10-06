# gui.py
# Tkinter GUI wrapper for building 'Fiber - Drop' GeoJSON from NAP + Service Locations.
# Uses modules.config.DATA_DIR for all internal I/O, per project rules.

import os
import shutil
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Always use modules.config.DATA_DIR for data I/O
from modules import config
from modules.drops.nap_to_sl_drop import build_drops_from_files


def get_downloads_dir() -> Path:
    # Cross-platform best-effort Downloads path
    home = Path.home()
    dl = home / "Downloads"
    return dl if dl.exists() else home


class DropsGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("NAP → Service Location Drops (Phase 1)")
        self.minsize(650, 260)

        # state
        self.nap_path_var = tk.StringVar()
        self.sl_path_var = tk.StringVar()
        self.out_dir_var = tk.StringVar(value=str(get_downloads_dir()))
        self.out_name_var = tk.StringVar(value="drops_phase1.geojson")
        self.output_chosen_by_user = False  # track if user explicitly set output dir

        # layout
        pad = {"padx": 10, "pady": 8}
        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True, **pad)

        # NAP file
        row = 0
        ttk.Label(frm, text="NAP GeoJSON:").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.nap_path_var).grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="Browse…", command=self.pick_nap).grid(row=row, column=2, sticky="e")

        # SL file
        row += 1
        ttk.Label(frm, text="Service Location GeoJSON:").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.sl_path_var).grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="Browse…", command=self.pick_sl).grid(row=row, column=2, sticky="e")

        # Output folder
        row += 1
        ttk.Label(frm, text="Output Folder:").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.out_dir_var).grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="Choose…", command=self.pick_out_dir).grid(row=row, column=2, sticky="e")

        # Output filename
        row += 1
        ttk.Label(frm, text="Output Filename:").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.out_name_var).grid(row=row, column=1, sticky="ew")
        ttk.Label(frm, text=".geojson").grid(row=row, column=2, sticky="w")

        # Run button
        row += 1
        run_btn = ttk.Button(frm, text="Build Drops", command=self.run_build)
        run_btn.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(12, 4))

        # Status
        row += 1
        self.status_var = tk.StringVar(value=f"DATA_DIR: {config.DATA_DIR}")
        ttk.Label(frm, textvariable=self.status_var, foreground="#555").grid(row=row, column=0, columnspan=3, sticky="w")

        # grid config
        frm.columnconfigure(1, weight=1)

        # if launched as .pyw, ensure GUI shows (this is mostly a no-op but explicit)
        try:
            if __file__.lower().endswith(".pyw"):
                self.deiconify()
        except Exception:
            pass

    # --- UI helpers ---

    def pick_nap(self):
        fp = filedialog.askopenfilename(
            title="Select NAP GeoJSON",
            filetypes=[("GeoJSON", "nap*.geojson"), ("All files", "*.*")],
            initialdir=self._initial_dir(),
        )
        if fp:
            self.nap_path_var.set(fp)
            self._maybe_set_output_dir(Path(fp).parent)

    def pick_sl(self):
        fp = filedialog.askopenfilename(
            title="Select Service Location GeoJSON",
            filetypes=[("GeoJSON", "service-location*.geojson"), ("All files", "*.*")],
            initialdir=self._initial_dir(),
        )
        if fp:
            self.sl_path_var.set(fp)
            self._maybe_set_output_dir(Path(fp).parent)

    def pick_out_dir(self):
        d = filedialog.askdirectory(
            title="Choose Output Folder",
            initialdir=self._initial_dir(),
        )
        if d:
            self.out_dir_var.set(d)
            self.output_chosen_by_user = True

    def _initial_dir(self) -> str:
        # Prefer last used folder among inputs; else Downloads
        for v in (self.nap_path_var.get(), self.sl_path_var.get()):
            if v:
                return str(Path(v).parent)
        return str(get_downloads_dir())

    def _maybe_set_output_dir(self, candidate: Path):
        # If user hasn't explicitly chosen an output directory, set it to input's folder
        if not self.output_chosen_by_user and candidate:
            self.out_dir_var.set(str(candidate))

    # --- core build flow ---

    def run_build(self):
        try:
            nap_src = Path(self.nap_path_var.get())
            sl_src  = Path(self.sl_path_var.get())

            if not nap_src.is_file() or not sl_src.is_file():
                messagebox.showerror("Missing Files", "Please select both NAP and Service Location GeoJSON files.")
                return

            out_dir = Path(self.out_dir_var.get().strip() or get_downloads_dir())
            out_dir.mkdir(parents=True, exist_ok=True)

            # Ensure DATA_DIR exists
            data_dir = Path(config.DATA_DIR)
            data_dir.mkdir(parents=True, exist_ok=True)

            # Copy inputs into DATA_DIR and call builder using basenames (per project rule)
            nap_dst = data_dir / nap_src.name
            sl_dst  = data_dir / sl_src.name
            shutil.copy2(nap_src, nap_dst)
            shutil.copy2(sl_src, sl_dst)

            # Run builder
            out_name = self.out_name_var.get().strip()
            if not out_name.lower().endswith(".geojson"):
                out_name += ".geojson"

            out_path = build_drops_from_files(
                nap_filename=nap_dst.name,
                service_locations_filename=sl_dst.name,
                out_filename=out_name,
            )

            # Copy the produced file from DATA_DIR to chosen output folder
            final_src = Path(out_path)
            final_dst = out_dir / final_src.name
            shutil.copy2(final_src, final_dst)

            self.status_var.set(f"Saved: {final_dst}")
            messagebox.showinfo("Success", f"Drops created:\n{final_dst}")

        except Exception as e:
            messagebox.showerror("Error", f"{type(e).__name__}: {e}")



if __name__ == "__main__":
    app = DropsGUI()
    app.mainloop()
