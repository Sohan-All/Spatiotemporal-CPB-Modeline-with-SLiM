import argparse
import os
import sys
import numpy as np
import pandas as pd

#!/usr/bin/env python3
"""
CreatePlots.py

Usage:
	python CreatePlots.py \
		--diversities diversities_2023.csv \
		--specifier specifier_matrix_2023.csv \
		[--out plot.png] [--cmap viridis] [--size 40]

This script reads a CSV of diversity values and a specifier CSV that gives
locations for those diversity values, then plots them on a plane using matplotlib.

Behavior:
- If the specifier CSV has exactly two numeric columns and the number of rows
	equals the number of diversity entries, it will plot a scatter using those
	two columns as x and y coordinates.
- Otherwise the script will try to interpret the specifier CSV as a matrix of
	integer indices that map into the diversity array (0- or 1-based). The mapped
	grid will be shown using imshow().
"""


import matplotlib.pyplot as plt


def read_diversities(path):
		df = pd.read_csv(path, header=0)
		# If there's a single column or a column named 'diversity' use that.
		if df.shape[1] == 1:
				values = df.iloc[:, 0].to_numpy(dtype=float)
		elif "diversity" in df.columns:
				values = df["diversity"].to_numpy(dtype=float)
		else:
				# fallback: use first numeric column
				num_cols = df.select_dtypes(include=[np.number]).columns
				if len(num_cols) == 0:
						raise ValueError("No numeric columns found in diversities CSV.")
				values = df[num_cols[0]].to_numpy(dtype=float)
		return values


def try_two_column_coords(spec_df, n_values):
		# Accept two numeric columns and number of rows matches values
		numeric_cols = spec_df.select_dtypes(include=[np.number]).columns
		if len(numeric_cols) >= 2 and spec_df.shape[0] == n_values:
				# specifier files are in the order: site, latitude, longitude
				# we want x=longitude, y=latitude for correct plotting
				x = spec_df[numeric_cols[1]].to_numpy(dtype=float)
				y = spec_df[numeric_cols[0]].to_numpy(dtype=float)
				return x, y
		return None


def try_matrix_mapping(spec_path, diversities):
		# Read as raw matrix (no header) to preserve integer matrix layout
		mat = pd.read_csv(spec_path, header=None).to_numpy()
		# If matrix is of float and matches diversity count when flattened, maybe it *is* the values
		flat = mat.flatten()
		if flat.size == diversities.size and np.issubdtype(flat.dtype, np.floating):
				# If values look identical to diversities (or plausible), use mat as grid of values
				return mat.astype(float), "direct_values"
		# If matrix contains integers in a suitable index range, map indices -> diversity values
		if np.issubdtype(mat.dtype, np.integer) or np.all(np.equal(np.mod(mat, 1), 0)):
				mat_int = mat.astype(int)
				n = diversities.size
				minv, maxv = mat_int.min(), mat_int.max()
				# 0-based indices
				if 0 <= minv and maxv < n:
						grid = np.vectorize(lambda i: diversities[int(i)])(mat_int)
						return grid, "index_0based"
				# 1-based indices
				if 1 <= minv and maxv <= n:
						grid = np.vectorize(lambda i: diversities[int(i) - 1])(mat_int)
						return grid, "index_1based"
		# If matrix shape flattened equals number of diversities, map by order
		if mat.size == diversities.size:
				grid = diversities.reshape(mat.shape)
				return grid, "reshape_in_order"
		# Otherwise cannot interpret as a matrix mapping
		return None, None


def plot_scatter(x, y, values, outpath, cmap, size, marker):
		fig, ax = plt.subplots(figsize=(8, 6))
		sc = ax.scatter(x, y, c=values, cmap=cmap, s=size, marker=marker, edgecolors='none')
		ax.set_xlabel("Longitude")
		ax.set_ylabel("Latitude")
		ax.set_title("Diversities")
		plt.colorbar(sc, ax=ax, label="Diversity (pi)")
		plt.tight_layout()
		fig.savefig(outpath, dpi=300)
		print(f"Saved scatter plot to {outpath}")


def plot_grid(grid, outpath, cmap, interpolation):
		fig, ax = plt.subplots(figsize=(8, 6))
		im = ax.imshow(grid, origin="lower", cmap=cmap, interpolation=interpolation)
		ax.set_title("Diversity grid")
		plt.colorbar(im, ax=ax, label="Diversity (pi)")
		plt.tight_layout()
		fig.savefig(outpath, dpi=300)
		print(f"Saved grid plot to {outpath}")


def main():
		parser = argparse.ArgumentParser(description="Plot diversities using a specifier matrix.")
		parser.add_argument("--diversities", default="diversities_2023.csv", help="CSV with diversity values")
		parser.add_argument("--specifier", default="specifier_matrix_2023.csv", help="CSV specifying locations or mapping")
		parser.add_argument("--out", default="diversity_plot.png", help="Output image path")
		parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
		parser.add_argument("--size", type=float, default=40.0, help="Point size for scatter")
		parser.add_argument("--marker", default="o", help="Marker for scatter")
		parser.add_argument("--interpolation", default="nearest", help="Interpolation for imshow")
		args = parser.parse_args()

		if not os.path.exists(args.diversities):
				print(f"Error: diversities file not found: {args.diversities}", file=sys.stderr)
				sys.exit(1)
		if not os.path.exists(args.specifier):
				print(f"Error: specifier file not found: {args.specifier}", file=sys.stderr)
				sys.exit(1)

		diversities = read_diversities(args.diversities)
		# Try reading specifier with headers first
		spec_df = pd.read_csv(args.specifier, header=0)
		coords = try_two_column_coords(spec_df, diversities.size)
		if coords is not None:
				x, y = coords
				plot_scatter(x, y, diversities, args.out, args.cmap, args.size, args.marker)
				return

		# Otherwise try matrix mapping heuristics
		grid, method = try_matrix_mapping(args.specifier, diversities)
		if grid is not None:
				plot_grid(grid, args.out, args.cmap, args.interpolation)
				return

		# As a last resort, if specifier CSV has same number of rows as diversities and at least 2 numeric cols, use first two numeric columns
		numeric_cols = spec_df.select_dtypes(include=[np.number]).columns
		if spec_df.shape[0] == diversities.size and len(numeric_cols) >= 2:
			# Ensure x is longitude and y is latitude
			x = spec_df[numeric_cols[1]].to_numpy(dtype=float)
			y = spec_df[numeric_cols[0]].to_numpy(dtype=float)
			plot_scatter(x, y, diversities, args.out, args.cmap, args.size, args.marker)
			return

		print("Could not interpret specifier CSV. Expected either:\n"
					"- a two-column x,y CSV with same number of rows as the diversities file, or\n"
					"- a matrix of integer indices mapping into the diversities vector (0- or 1-based), or\n"
					"- a matrix whose shape matches the desired grid and whose cells are the diversity values.",
					file=sys.stderr)
		sys.exit(2)


if __name__ == "__main__":
		main()