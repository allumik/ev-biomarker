import pandas as pd
import numpy as np
import seaborn as sns
import altair as alt
import matplotlib.pyplot as plt
from typing import Tuple
from scipy.stats import spearmanr
from scipy.stats import f_oneway
from scipy.stats import chi2
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.decomposition import PCA

import dendro_barplotter


# %% Statistical testing functoins
def convert_to_long(fracs, pheno) -> pd.DataFrame:
  # convert to long format
  return (
  fracs
    .reset_index(names="index")
    .melt(id_vars="index", var_name="celltype", value_name="ratio")
    .set_index("index")
    .merge(pheno, left_index=True, right_index=True)
  )

# Perform one-way ANOVA and Tukey's test as a function for grouping
def pairwise_tukey_df(dat) -> pd.DataFrame:
  tukey_res = pairwise_tukeyhsd(
      dat["ratio"],
      dat["cyclephase"],
    )._results_table
  return pd.DataFrame(
    data=tukey_res.data[1:],
    columns=tukey_res.data[0]
  )

def tukey_table(fracs_long) -> pd.DataFrame:
  tukey_res = {
    celltype:
      pairwise_tukey_df(dat)
    for celltype, dat in fracs_long.groupby("celltype")
  }
  return pd.concat(
    tukey_res, keys=tukey_res.keys(), names=["celltype"]
    ).reset_index(level='celltype')

def anova_table(fracs_long) -> pd.DataFrame:
  anova_res = {
    celltype:
      f_oneway(*[dat.query("cyclephase == @cycle")["ratio"] for cycle in dat.cyclephase.unique()])
    for celltype, dat in fracs_long.groupby("celltype")
  }
  return pd.DataFrame(anova_res).T.set_axis(["f-stat", "p-val"], axis=1)

def cor_table(fracs_long) -> pd.DataFrame:
  # convert the cyclephase to a dummy var
  fracs_long["dummy_phase"] = pd.Categorical(
    fracs_long.cyclephase,
    categories=["pre", "pro", "rec", "post"],
    ordered=True
    ).codes
  cor_res = {
    celltype:
      spearmanr(dat["dummy_phase"], dat["ratio"])
    for celltype, dat in fracs_long.groupby("celltype")
  }
  return pd.DataFrame(cor_res).T.set_axis(["corr", "p-val"], axis=1)

def mahalanobis_distances(x, refs):
  cov_val = np.cov(refs)
  inv_cov = 1 / cov_val
  diff = x - np.mean(refs)
  mahal = np.sqrt(diff**2 * inv_cov)
  return mahal

# Function to calculate geometric mean, ignoring NaN and non-positive values
def geometric_mean(row):
    valid_values = row[(row > 0) & (abs(row) < np.inf)].dropna()
    if len(valid_values) > 0:
        return np.exp(np.log(valid_values).mean())
    else:
        return np.nan



# %% Functions for plotting
def dendro_barplot(
  fractions: pd.DataFrame,
  phenotype_df: pd.DataFrame,
  legend_phase_bbox: tuple = (1, 0.3),
  legend_class_bbox: tuple = (1, 0.3),
  legend_phase_ncol: int = 1,
  legend_class_ncol: int = 1,
  phenotype_df_col: str = "cyclephase",
  add_dendrogram: bool = True
) -> None:
  """Generates a dendrogram-barplot figure to visualize cell fraction predictions.

  This function creates a combined plot consisting of a dendrogram (optionally)
  and a stacked bar plot. The dendrogram visualizes the hierarchical clustering
  of samples based on their cell fractions. The bar plot shows the predicted
  cell fractions for each sample, with bars stacked according to different
  cell types/states. The rows of the bar plot are colored according to
  phenotype data.

  Args:
      fractions: DataFrame of cell fraction predictions. Rows are samples,
          columns are cell types/states. Values are the predicted fractions.
      phenotype_df: DataFrame of phenotype data. Must have the same index as `fractions`.
      legend_phase_bbox: Bounding box coordinates (x, y) for the phase legend.
          This controls the legend's position.
      legend_class_bbox: Bounding box coordinates (x, y) for the class legend.
      legend_phase_ncol: Number of columns for the phase legend.
      legend_class_ncol: Number of columns for the class legend.
      phenotype_df_col: Column name in `phenotype_df` to use for row
          coloring/labels (e.g., 'cyclephase'). This column determines
          the categories used for coloring the rows of the bar plot.
      add_dendrogram: Whether to include the dendrogram in the plot.
  """

  frac_pred_index = fractions.index

  # 1. Create colormap and assign row colors
  lut_row = dendro_barplotter.create_row_colormap(phenotype_df, phenotype_df_col)
  row_colors = dendro_barplotter.assign_row_colors(phenotype_df, phenotype_df_col, frac_pred_index, lut_row)

  # 2. Setup figure and axes
  fig, (ax_dendro, ax_bar) = dendro_barplotter.setup_figure_axes(add_dendrogram)

  dend_obj = None  # Initialize dend_obj with type hint

  # 3. Plot dendrogram and style axes (conditional)
  if add_dendrogram and ax_dendro is not None:  # Check if ax_dendro is valid
      dend_obj = dendro_barplotter.plot_dendrogram(ax_dendro, fractions, frac_pred_index, phenotype_df, phenotype_df_col)
      dendro_barplotter.style_dendrogram_axes(ax_dendro)

  # 4. Plot barplot and style axes
  dendro_barplotter.plot_barplot(ax_bar, fractions, dend_obj, row_colors)  # Pass dend_obj to barplot
  dendro_barplotter.style_barplot_axes(ax_bar)

  # 5. Create legends
  dendro_barplotter.create_phase_legend(ax_bar, lut_row, phenotype_df_col, legend_phase_bbox, legend_phase_ncol)
  dendro_barplotter.create_class_legend(ax_bar, legend_class_bbox, legend_class_ncol)

  plt.show()



def biplot_fractions_altair(
  fractions_df: pd.DataFrame,
  phenotype_df: pd.DataFrame,
  legend_title: str = "Cycle Phase",
  color_field: str = None,  # Added for flexibility
  style_field: str = None   # Added for flexibility
) -> Tuple[alt.Chart, PCA]:
  """Generates an Altair biplot of Principal Component Analysis (PCA) and its loadings.

  Performs PCA, creates an interactive scatter plot, and overlays loadings. Returns the Altair chart
  and the fitted PCA object.

  Args:
      fractions_df: DataFrame where rows are samples, columns are features (fractions).
      phenotype_df: DataFrame with phenotype data. Must have the same index as `fractions_df`.
      legend_title: Title for the legend.
      color_field: Column in the merged DataFrame to use for color encoding.
      style_field: Column in the merged DataFrame to use for shape encoding.

  Returns:
      Tuple[alt.Chart, PCA]: A tuple containing:
          - alt.Chart: The Altair chart object representing the biplot.
          - PCA: The fitted scikit-learn PCA object.
  """

  # PCA calculation
  pca = PCA(n_components=2)
  pca_df = (
    pd.DataFrame(data=pca.fit_transform(fractions_df), columns=["PC1", "PC2"])
    .set_index(fractions_df.index)
    .merge(phenotype_df, left_index=True, right_index=True)
  )
  expl_var = pca.explained_variance_ratio_
  loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
  loadings_df = pd.DataFrame(
    loadings, index=fractions_df.columns, columns=["PC1", "PC2"]
  ).reset_index()  # Reset index for Altair

  # --- Create the base scatter plot ---
  base = alt.Chart(pca_df).encode(
    x=alt.X("PC1:Q", axis=alt.Axis(title=f"PC1 ({expl_var[0]:.2%})")),
    y=alt.Y("PC2:Q", axis=alt.Axis(title=f"PC2 ({expl_var[1]:.2%})")),
  )

  # Add color and shape encoding if specified
  if color_field:
    base = base.encode(color=alt.Color(f"{color_field}:N", title=legend_title))
  if style_field:
    base = base.encode(shape=alt.Shape(f"{style_field}:N"))
  
  scatter = base.mark_point(size=60)

  # --- Create the loadings plot ---
  loadings_chart = (
    alt.Chart(pd.concat([
      loadings_df, 
      pd.DataFrame({"index": loadings_df["index"], "PC1":0.0, "PC2": 0})
      ], ignore_index=True))
    .mark_line(color="red", opacity=0.5)  # Use mark_rule instead of mark_line
    .encode(
      x='PC1:Q',
      y='PC2:Q',
      detail="index"
    )
  )

  # --- Create text labels for loadings (with filtering) ---
  loadings_text = (
    alt.Chart(loadings_df)
      .mark_text(align='left', dx=5, dy=-5, color="black", fontSize=10)
      .encode(
        x="PC1:Q",
        y="PC2:Q",
        text='index:N'
      )
      .transform_filter(
      (alt.datum.PC1 > 0.01) | (alt.datum.PC1 < -0.01) | (alt.datum.PC2 > 0.01) | (alt.datum.PC2 < -0.01)
      )
  )

  # --- Combine the plots ---
  chart = (scatter + loadings_chart + loadings_text).properties(
    title="PCA Biplot",
    width=600,  # Adjust as needed
    height=400   # Adjust as needed
  )

  return chart, pca


def biplot_fractions(fractions_df: pd.DataFrame, phenotype_df: pd.DataFrame, legend_title: str = "Cycle Phase", **kwargs) -> Tuple[plt.Figure, PCA]:
  """Generates a biplot of Principal Component Analysis (PCA) and its loadings.

  This function performs PCA on the input `fractions_df` (DataFrame representing
  proportions or fractions of different features), and visualizes the results
  as a scatter plot (biplot).  It overlays the loadings of the original features
  onto the PCA plot, showing their contribution to the principal components.

  Args:
      fractions_df: A pandas DataFrame where rows represent samples and columns represent features.
          The values in the DataFrame should be fractions or proportions (not raw counts).
      phenotype_df: A pandas DataFrame containing phenotype data. It must have the same index as `fractions_df`.
          This DataFrame is used to color or otherwise style the points in the scatter plot,
          based on the provided keyword arguments.
      legend_title: The title for the legend in the plot. Defaults to "Cycle Phase".
      **kwargs:  Additional keyword arguments passed to `seaborn.scatterplot`.  This allows
          customization of the scatter plot, such as coloring points by a column in `phenotype_df`
          (e.g., `hue=phenotype_df["Cycle Phase"]`).
  """
  # PCA calculation
  pca = PCA(n_components=2)
  pca_df = (
    pd.DataFrame(data=pca.fit_transform(fractions_df), columns=["PC1", "PC2"])
    .set_index(fractions_df.index)
    .merge(phenotype_df, left_index=True, right_index=True)
  )
  expl_var = pca.explained_variance_ratio_
  loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

  # Visualization
  fig = plt.figure(figsize=(8, 6))
  sns.scatterplot(x='PC1', y='PC2', data=pca_df, **kwargs)

  # Plot loadings
  for i, feature in enumerate(fractions_df.columns):  # Corrected: Use frac_obj.columns
      plt.arrow(0, 0, loadings[i, 0], loadings[i, 1], color='r', alpha=0.5, head_width=0.005, linewidth=.6)
      # Display text for significant loadings
      if abs(loadings[i, 0]) > 0.01 or abs(loadings[i, 1]) > 0.01:  #Corrected to use abs
          plt.text(loadings[i, 0] * 1.15, loadings[i, 1] * 1.4, feature, color='black', ha='center', va='center', size=8)

  plt.title('PCA Result')
  plt.xlabel(f'PC1 ({expl_var[0]:.2%})')
  plt.ylabel(f'PC2 ({expl_var[1]:.2%})')
  plt.legend(title=legend_title, loc="upper left", bbox_to_anchor=(1.05, 1))
  plt.grid(False)
  # Major grid only
  plt.axhline(0, color="lightgray", linewidth=1)
  plt.axvline(0, color="lightgray", linewidth=1)

  return fig, pca



def peruvian_transform(
  fractions: pd.DataFrame,
  phenotype: pd.DataFrame,
  general_cells: pd.DataFrame = None,
  grouping_ids: list = None
  ) -> pd.DataFrame:
  """Transforms and reshapes data for stacked area plots of cell type fractions.

  Prepares a DataFrame for visualization (e.g., with Altair) by:
    1. Optionally grouping and averaging cell type fractions by specified columns.
    2. Reshaping the data from wide to long format (stacking cell types).
    3. Optionally merging cell lineage information based on 'celltype'.

  Args:
      fractions: DataFrame of cell type fraction predictions.  Rows are samples,
          columns are cell types.  Should contain columns to join with `phenotype`.
      phenotype: DataFrame of phenotype data. Must contain the columns specified
          in `grouping_ids`, and index or columns to join with `fractions`.
      general_cells: DataFrame mapping cell types to lineages.
          Requires 'celltype' and 'lineage' columns. If None, no lineage information is added.
      grouping_ids: List of column names in `phenotype` to group by.
          Defaults to ["cyclephase"].  If an empty list is provided (`[]`), no grouping will occur
          and all samples will be present for the plot.

  Returns:
      DataFrame: Transformed DataFrame in long format, ready for plotting.
          Columns include:
          - The columns listed in `grouping_ids`.
          - 'celltype' (the cell type).
          - 'fractions' (the fraction of that cell type).
          - 'lineage' (if `general_cells` is provided).
  """
  # Group and average the samples
  if grouping_ids is not None:
    fractions = fractions.join(phenotype[grouping_ids]).groupby(grouping_ids).mean()

  # Data Reshaping for Altair: Convert to long format
  fractions.columns.names = ["celltype"]
  fractions = (
    fractions
    .stack()
    .reset_index()
    .rename(columns={0: "fractions"})
    .reset_index(drop=True)
  )
  # Add lineage information for the celltype
  if general_cells is not None:  # More concise conditional merge
      return fractions.merge(general_cells, on="celltype")
  else:
      return fractions # Add lineage information for the celltype



def peruvian_sands(
  fractions_df: pd.DataFrame,
  x_axis_col: str = "cyclephase",
  x_order: list = ["pro", "pre", "rec", "post"],
  x_title: str = None,
  global_color_scale: alt.Scale = None,
  lineage_val: str = None,
  grouping_val: str = None,
  is_first_dataset_in_row: bool = True,
  is_first_row: bool = True,
  legend_col_n: int = 1,
  dims: tuple = (200, 150)
) -> alt.Chart:
  """Generates a stacked area chart of cell type fractions across cycle phases.

  Visualizes cell type fractions over cycle phases using Altair-Lite.
  Designed for faceted/concatenated layouts, conditionally displaying
  axis titles, legends, and plot titles based on layout position.

  Args:
      fractions_df: DataFrame with 'cyclephase', 'fractions', and 'celltype' columns.
      cyclephase_order: Order of cycle phases for x-axis.
      lineage_val: Lineage for y-axis title (if first in row).
      grouping_val: Dataset for plot title (if first row).
      is_first_dataset_in_row: If True, display y-axis title and legend.
      is_first_row: If True, display plot title.
      dims: (width, height) of the plot in pixels.

  Returns:
      An Altair-Lite Chart object.
  """
  if global_color_scale is None: global_color_scale = alt.Scale(domain=fractions_df.celltype.unique().tolist())
  return alt.Chart(fractions_df).mark_area().encode(
    x=alt.X(
      f'{x_axis_col}:O',
      axis=alt.Axis(labelAngle=-45, grid=False, domain=False, title=x_title),
      sort=x_order
      ),
    y=alt.Y(
      'fractions:Q',
      axis=alt.Axis(format='%', grid=False, title=lineage_val) if is_first_dataset_in_row else None
      ),
    color=alt.Color(
      'celltype:N', # Show legend only for the first chart in the row
      scale=global_color_scale,
      legend=alt.Legend(
        title=None,
        values=fractions_df.celltype.unique().tolist(), # Assuming current_row_celltypes is defined in the outer scope or needs to be passed if row-specific
        columns=legend_col_n
        ) if is_first_dataset_in_row else None
      ),
  ).properties( # Title only for first plot, or dataset name as subtitle
    width=dims[0],
    height=dims[1],
    title=f'{grouping_val}' if is_first_row or not None else ""
  )



def peruvian_grouped(
  fractions_df: pd.DataFrame,
  cyclephase_order: list = ["pro", "pre", "rec", "post"],
  grouping_id: str = "dataset",
  dims: tuple = (200, 150)
) -> alt.Chart:
  """Generates a faceted stacked area chart of cell type fractions across cycle phases, faceted by lineage and dataset.

  This function takes cell type fraction predictions, phenotype data, and cell lineage information to create a
  faceted stacked area chart using Altair-Lite. The chart visualizes cell type fractions across different
  cycle phases, with facets arranged by cell lineage (rows) and grouping ID (columns, typically dataset).
  Each row (lineage) has a shared y-axis scale, and the entire grid shares an x-axis scale. A global color
  scheme is used for cell types, and legends are displayed only at the beginning of each row, showing
  cell types present in that lineage.

  Args:
    fractions_df (pd.DataFrame): DataFrame containing cell type fraction predictions.
      Must have columns that will be used to join with `comb_uf_pheno` based on 'sample_id',
      and columns representing cell types (whose names will become the columns of the 'fractions' table
      after grouping and averaging).
    comb_uf_pheno (pd.DataFrame): DataFrame containing phenotype data, including 'cyclephase', and the
      column specified by `grouping_id` (e.g., 'dataset'). Must have a 'sample_id' column to join with `frac_pred`.
    general_cells (pd.DataFrame): DataFrame containing cell type to lineage mapping.
      Must have 'celltype' and 'lineage' columns.
    grouping_id (str, optional): The column name in `comb_uf_pheno` to use for grouping datasets
      horizontally in the faceted plot. Defaults to "dataset".

  Returns:
    Chart: An Altair-Lite Chart object representing the faceted stacked area chart.
  """
  row_charts = [] # List to hold charts for each row (lineage)
  is_first_row = True
  global_color_scale = alt.Scale(domain=fractions_df.celltype.unique().tolist())

  for lineage_val in fractions_df.lineage.unique():
    dataset_charts_row = [] # List to hold charts for each dataset in the current lineage row
    is_first_dataset_in_row = True # Flag to indicate if it's the first plot in the row (for legend)

    for grouping_val in fractions_df[grouping_id].unique():
      # Create a chart for this dataset and lineage
      current_chart = peruvian_sands(
        fractions_df.query(f"lineage == '{lineage_val}' & {grouping_id} == '{grouping_val}'"),
        x_order=cyclephase_order,
        lineage_val=lineage_val,
        grouping_val=grouping_val,
        is_first_dataset_in_row=is_first_dataset_in_row,
        is_first_row=is_first_row,
        global_color_scale=global_color_scale, # Pass global color scale
        dims=dims
      )
      dataset_charts_row.append(current_chart)
      is_first_dataset_in_row = False
    is_first_row = False

    # Horizontally concatenate charts for the current lineage row
    row_chart = alt.hconcat(*dataset_charts_row, spacing=0).resolve_scale(y='shared')
    row_charts.append(row_chart) # Add the row chart to the list of rows

  # Vertically concatenate all row charts to form the final grid
  return alt.vconcat(*row_charts, spacing=0).resolve_scale( # Adjust vertical spacing as needed
    x='shared',
    y='independent',
    color='independent' # Keep color scale independent - legends are now handled manually
  ).properties(
    config=alt.Config(
      view=alt.ViewConfig(
        strokeWidth=0 # Set strokeWidth to 0 to remove the border line
      )
    )
  )