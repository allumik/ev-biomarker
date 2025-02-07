import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import f_oneway
from scipy.stats import chi2
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


def merge_and_group_data(fractions, pheno_dat, pheno_dat_cols):
  """
  Merges fractions and phenotype data and groups by specified columns.

  Args:
      fractions (pd.DataFrame): DataFrame of cell type fractions, index is sample names.
    pheno_dat (pd.DataFrame): Phenotype data DataFrame, index is sample names.
    pheno_dat_cols (list): List of column names from pheno_dat to group by.

  Returns:
    pd.DataFrame: DataFrame of average fractions grouped by pheno_dat_cols.
  """
  merged_data = pd.merge(fractions, pheno_dat[pheno_dat_cols], left_index=True, right_index=True)
  grouped_avg_fractions = merged_data.groupby(pheno_dat_cols).mean().reset_index()
  return grouped_avg_fractions

def create_cyclephase_colormap(pheno_dat, pheno_dat_cols):
  """
  Creates a colormap for cycle phases if 'cyclephase' is in pheno_dat_cols.

  Args:
    pheno_dat (pd.DataFrame): Phenotype data DataFrame.
    pheno_dat_cols (list): List of column names used for grouping, may include 'cyclephase'.

  Returns:
    dict: Colormap for cycle phases, empty dict if 'cyclephase' is not in pheno_dat_cols or pheno_dat.columns.
  """
  if 'cyclephase' in pheno_dat_cols and 'cyclephase' in pheno_dat.columns:
    return dict(zip(
      pheno_dat['cyclephase'].unique(),
      sns.color_palette("Set2", pheno_dat['cyclephase'].unique().size).as_hex()
    ))
  return {}

def create_barplot(grouped_avg_fractions, pheno_dat_cols, lut_row, ax_bar): # Modified to accept ax_bar
  """
  Creates the barplot using seaborn on a provided axes.

  Args:
      grouped_avg_fractions (pd.DataFrame): DataFrame of average fractions.
      pheno_dat_cols (list): List of column names used for grouping.
      lut_row (dict): Colormap for cycle phases.
      ax_bar (matplotlib.axes._axes.Axes): Axes object to plot on. # Added ax_bar docstring

  Returns:
      matplotlib.axes._axes.Axes: Axes object with the created barplot. # Modified return docstring
  """

  plot_data_melted = pd.melt(grouped_avg_fractions,
                      id_vars=pheno_dat_cols,
                      var_name='CellType',
                      value_name='Average Fraction')
  plot_data_melted['GroupLabel'] = plot_data_melted[pheno_dat_cols].apply(lambda row: '_'.join(row.astype(str)), axis=1)

  sns.barplot(x='GroupLabel', y='Average Fraction', hue='CellType', data=plot_data_melted, ax=ax_bar) # Removed stacked=True, plotting on provided ax_bar

  group_label_string = ', '.join([col.capitalize() for col in pheno_dat_cols])
  ax_bar.set_xlabel(group_label_string)
  ax_bar.set_ylabel('Average Cell Fraction')
  ax_bar.set_title('Average Cell Fraction by ' + group_label_string)

  minimal_barplot_theme(ax_bar) # Apply minimal theme
  if 'cyclephase' in pheno_dat_cols and lut_row:
    color_xtick_labels_by_cyclephase(ax_bar, lut_row) # Color x tick labels if cyclephase is used
  else:
    ax_bar.set_xticklabels(ax_bar.get_xticklabels(), rotation=45, ha='right') # Still rotate if no cyclephase coloring

  return ax_bar # Return the modified axes


def minimal_barplot_theme(ax):
  """
  Applies a minimal theme to the barplot axes (removes spines and grid).

  Args:
    ax (matplotlib.axes._axes.Axes): Barplot axes object.
  """
  for spine in ax.spines.keys():
    ax.spines[spine].set_visible(False)
  ax.grid(False)

def color_xtick_labels_by_cyclephase(ax_bar, lut_row):
  """
  Colors x-tick labels of the barplot based on cycle phase information.

  Args:
      ax_bar (matplotlib.axes._axes.Axes): Barplot axes object.
      lut_row (dict): Colormap for cycle phases.
  """
  xtick_labels = ax_bar.get_xticklabels()
  for label in xtick_labels:
    group_label_parts = label.get_text().split('_')
    for part in group_label_parts:
      if part in lut_row:
        label.set_color(lut_row[part])
        break
  ax_bar.set_xticklabels(xtick_labels, rotation=45, ha='right')


def create_cyclephase_legend(ax_bar, lut_row, pheno_dat_cols, legend_bbox, legend_ncol):
  """
  Creates the cycle phase legend for the barplot if 'cyclephase' is used for grouping.

  Args:
  ax_bar (matplotlib.axes._axes.Axes): Barplot axes object.
  lut_row (dict): Colormap for cycle phases.
  pheno_dat_cols (list): List of column names used for grouping.
  legend_bbox (tuple): Bounding box for the legend.
  legend_ncol (int): Number of columns in the legend.

  Returns:
  matplotlib.legend.Legend or None: Cycle phase legend object, None if no legend is created.
  """
  if 'cyclephase' in pheno_dat_cols and lut_row:
    legend_elements_phases = [
    Patch(facecolor=color, edgecolor='black', label=phase) for phase, color in lut_row.items()
    ]
    legend_phases = ax_bar.legend(
      handles=legend_elements_phases,
      title='Cyclephase',
      bbox_to_anchor=legend_bbox,
      loc='upper left',
      ncol=legend_ncol,
      frameon=False
    )
    return legend_phases
  return None # Return None if no legend is created

def create_celltype_legend(ax_bar, legend_bbox, legend_ncol):
  """
  Creates the cell type legend for the barplot.

  Args:
  ax_bar (matplotlib.axes._axes.Axes): Barplot axes object.
  legend_bbox (tuple): Bounding box for the legend.
  legend_ncol (int): Number of columns in the legend.

  Returns:
  matplotlib.legend.Legend: Cell type legend object.
  """
  handles, labels = ax_bar.get_legend_handles_labels()
  legend_classes = ax_bar.legend(
    handles=handles,
    labels=labels,
    title='Cell Types',
    bbox_to_anchor=legend_bbox,
    loc='lower left',
    ncol=legend_ncol,
    frameon=False
  )
  return legend_classes


# --- Main Function ---
def grouped_barplot_multiple_pheno_cols_average_modular(
    fractions,
    pheno_dat,
    pheno_dat_cols=["cyclephase"],
    legend_phase_bbox=(1, .3),
    legend_class_bbox=(1, .3),
    legend_phase_ncol=1,
    legend_class_ncol=1,
    figsize=(10, 6)
    ):
  """
  Generates a grouped average barplot of cell type fractions.

  The barplot groups samples based on the columns specified in `pheno_dat_cols`
  and displays the average cell type fractions for each group.

  If 'cyclephase' is included in `pheno_dat_cols`, the function will:
  1. Color the x-tick labels of the barplot according to the cycle phase.
  2. Include a legend for cycle phases, mapping colors to phase names.

  Args:
    fractions (pd.DataFrame): DataFrame of cell type fractions.
        Index should be sample names, and columns should be cell type names.
    pheno_dat (pd.DataFrame): Phenotype data DataFrame.
        Index should be sample names, and must contain columns listed in `pheno_dat_cols`.
    pheno_dat_cols (list): List of column names from `pheno_dat` to group samples by.
        Defaults to `["cyclephase"]`. Can be a list of multiple columns for multi-level grouping.
    legend_phase_bbox (tuple): Bounding box coordinates (x, y) for the cycle phase legend
        in axes coordinates. Defaults to `(1, .3)` (placed outside the plot to the right).
    legend_class_bbox (tuple): Bounding box coordinates (x, y) for the cell type legend.
        Defaults to `(1, .3)` (placed outside the plot to the right).
    legend_phase_ncol (int): Number of columns in the cycle phase legend. Defaults to 1.
    legend_class_ncol (int): Number of columns in the cell type legend. Defaults to 1.
    figsize (tuple): Figure size (width, height) in inches. Defaults to `(10, 6)`.

  Returns:
      None: Displays the plot using `plt.show()`.
  """

  grouped_avg_fractions = merge_and_group_data(fractions, pheno_dat, pheno_dat_cols)
  lut_row = create_cyclephase_colormap(pheno_dat, pheno_dat_cols)
  fig, ax_bar = create_barplot(grouped_avg_fractions, pheno_dat_cols, lut_row, figsize)

  legend_phases = create_cyclephase_legend(ax_bar, lut_row, pheno_dat_cols, legend_phase_bbox, legend_phase_ncol)
  legend_classes = create_celltype_legend(ax_bar, legend_class_bbox, legend_class_ncol)

  if legend_phases:
    ax_bar.add_artist(legend_phases)

  plt.tight_layout()
  plt.show()



# %% Plot a barplot with a dendrogram on top of it
# TODO: refactor it to use the previous functoins


def dendro_barplot(
  fractions,
  pheno_dat,
  pheno_dat_cols=["cyclephase"],
  legend_phase_bbox=(1, .3),
  legend_class_bbox=(1, .3),
  legend_phase_ncol=1,
  legend_class_ncol=1,
  figsize=(12, 10)
  ):
  """
  Generates a dendrogram combined with a grouped average barplot of cell type fractions.

  The barplot groups samples based on the columns specified in `pheno_dat_cols`
  and displays the average cell type fractions for each group. If 'cyclephase'
  is in `pheno_dat_cols`, the x-tick labels are colored by cycle phase, and a
  cycle phase legend is added.

  Args:
      fractions (pd.DataFrame): DataFrame of cell type fractions, index is sample names,
        columns are cell type names.
      pheno_dat (pd.DataFrame): Phenotype data DataFrame, index is sample names,
        contains columns specified in `pheno_dat_cols`.
      pheno_dat_cols (list): List of column names from pheno_dat to group samples by in the barplot.
        Defaults to ["cyclephase"].
      legend_phase_bbox (tuple): Bounding box coordinates for the cycle phase legend.
      legend_class_bbox (tuple): Bounding box coordinates for the cell type legend.
      legend_phase_ncol (int): Number of columns in the cycle phase legend.
      legend_class_ncol (int): Number of columns in the cell type legend.
      figsize (tuple): Figure size for the entire plot.

  Returns:
      None: Displays the plot using plt.show().
  """
  frac_pred_index = fractions.index

  # 1. Create a colormap for rows (cycle phases) - using modular function
  lut_row = create_cyclephase_colormap(pheno_dat, pheno_dat_cols)
  row_colors = None # Not using row_colors directly for bar coloring anymore in this version

  # 2. Set up figure with dendrogram and barplot axes
  fig = plt.figure(figsize=figsize)
  gs = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)
  ax_dendro = fig.add_subplot(gs[0])
  ax_bar = fig.add_subplot(gs[1])

  # 3. Dendrogram part (same as before)
  linkage_mat = linkage(fractions.loc[frac_pred_index], method="complete", optimal_ordering=True)
  dendrogram(
      linkage_mat,
      labels=pheno_dat[pheno_dat_cols[0]][frac_pred_index].index if pheno_dat_cols and pheno_dat_cols[0] in pheno_dat else pheno_dat[frac_pred_index].index, # Use first pheno_col or index as labels
      color_threshold=0,
      above_threshold_color="#000000",
      ax=ax_dendro
      )
  ax_dendro.set_title('TAPE deconvolution results')
  # Minimal theme for dendrogram
  for spine in ax_dendro.spines.keys(): ax_dendro.spines[spine].set_visible(False)
  ax_dendro.tick_params(left=False, bottom=False)
  ax_dendro.set_xticks([])
  ax_dendro.set_yticks([])

  # 4. Grouped average barplot part - using modular functions
  grouped_avg_fractions = merge_and_group_data(fractions, pheno_dat, pheno_dat_cols)
  create_barplot(grouped_avg_fractions, pheno_dat_cols, lut_row, ax_bar) # Create plot using modular function, pass ax_bar

  # 5. Legends - using modular legend functions
  legend_phases = create_cyclephase_legend(ax_bar, lut_row, pheno_dat_cols, legend_phase_bbox, legend_phase_ncol)
  legend_classes = create_celltype_legend(ax_bar, legend_class_bbox, legend_class_ncol)

  if legend_phases:
    ax_bar.add_artist(legend_phases)

  plt.tight_layout() # Use tight layout to adjust for both subplots
  plt.show()


# %% To plot the biplot
def biplot_fractions(frac_obj, pheno, legend_title="Cycle Phase", kwargs={}):
  ## don't scale as the counts are props anyways
  pca = PCA(n_components=2)
  pca_df = (
  pd.DataFrame(data=pca.fit_transform(frac_obj), columns=["PC1", "PC2"])
  .set_index(frac_obj.index)
  .merge(pheno, left_index=True, right_index=True)
  )
  expl_var = pca.explained_variance_ratio_
  loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

  # Visualize the PCA results
  plt.figure(figsize=(8, 6))
  sns.scatterplot(x='PC1', y='PC2', data=pca_df, **kwargs)

  # Plot the loadings as a biplot
  for i, feature in enumerate(frac_obj):
    plt.arrow(0, 0, loadings[i, 0], loadings[i, 1], color='r', alpha=0.5, head_width=0.005, linewidth=.6)
    # write only signf loadings
    if loadings[i, 0] > 0.01 or loadings[i, 1] > 0.01:
      plt.text(loadings[i, 0] * 1.15, loadings[i, 1] * 1.4, feature, color='black', ha='center', va='center', size=8)

  plt.title('PCA Result')
  plt.xlabel(f'PC1 ({expl_var[0]:.2%})')
  plt.ylabel(f'PC2 ({expl_var[1]:.2%})')
  plt.legend(title=legend_title, loc="upper left", bbox_to_anchor=(1.05, 1))
  plt.grid(False)
  # Major grid only
  plt.axhline(0, color="lightgray", linewidth=1)
  plt.axvline(0, color="lightgray", linewidth=1)
  plt.show()



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
