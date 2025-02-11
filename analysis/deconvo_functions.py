import pandas as pd
import numpy as np
import seaborn as sns
import altair as alt
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import f_oneway
from scipy.stats import chi2
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.decomposition import PCA

import dendro_barplotter



def dendro_barplot(
    fractions,
    pheno_dat,
    legend_phase_bbox=(1, .3),
    legend_class_bbox=(1, .3),
    legend_phase_ncol=1,
    legend_class_ncol=1,
    pheno_dat_col="cyclephase",
    add_dendro=True
    ):
  """
  Generates a dendrogram-barplot figure to visualize cell fraction predictions.

  Args:
      fractions (pd.DataFrame): DataFrame of cell fraction predictions.
      pheno_dat (pd.DataFrame): DataFrame of phenotype data.
      legend_phase_bbox (tuple): Bounding box for the phase legend.
      legend_class_bbox (tuple): Bounding box for the class legend.
      legend_phase_ncol (int): Number of columns in the phase legend.
      legend_class_ncol (int): Number of columns in the class legend.
      pheno_dat_col (str): Column in pheno_dat to use for row coloring/labels.
      add_dendro (bool): Whether to include the dendrogram in the plot (default: True).
  """
  frac_pred_index = fractions.index

  # 1. Create colormap and assign row colors
  lut_row = dendro_barplotter.create_row_colormap(pheno_dat, pheno_dat_col)
  row_colors = dendro_barplotter.assign_row_colors(pheno_dat, pheno_dat_col, frac_pred_index, lut_row)

  # 2. Setup figure and axes
  fig, (ax_dendro, ax_bar) = dendro_barplotter.setup_figure_axes(add_dendro)

  dend_obj = None # Initialize dend_obj to None

  # 3. Plot dendrogram and style axes (conditional)
  if add_dendro and ax_dendro is not None: # Check if ax_dendro is valid
    dend_obj = dendro_barplotter.plot_dendrogram(ax_dendro, fractions, frac_pred_index, pheno_dat, pheno_dat_col)
    dendro_barplotter.style_dendrogram_axes(ax_dendro)

  # 4. Plot barplot and style axes
  dendro_barplotter.plot_barplot(ax_bar, fractions, dend_obj, row_colors) # Pass dend_obj to barplot
  dendro_barplotter.style_barplot_axes(ax_bar)

  # 5. Create legends
  dendro_barplotter.create_phase_legend(ax_bar, lut_row, pheno_dat_col, legend_phase_bbox, legend_phase_ncol)
  dendro_barplotter.create_class_legend(ax_bar, legend_class_bbox, legend_class_ncol)

  plt.show()



# TODO: rewrite the barplot thing to use sankey plots for the proprtions.
def sankey_props(fractions, pheno_dat, grouping_col):
  # 1. Data Reshaping for Altair: Convert to long format
  df_long = fractions.stack().reset_index()
  df_long.columns = ['Sample', 'Cell_Type', 'Fraction'] # Meaningful column names

  # 2. Create the Stacked Area Chart with Altair-Lite
  chart = alt.Chart(df_long).mark_area().encode(
      x=alt.X('Sample:N',  # Nominal (categorical) data for samples
              axis=alt.Axis(labelAngle=-45)), # Rotate x-axis labels for readability
      y=alt.Y('Fraction:Q',  # Quantitative data for fraction
              axis=alt.Axis(format='%')), # Format y-axis as percentages
      color=alt.Color('Cell_Type:N',  # Nominal data for cell types, for stacking and color differentiation
                      legend=alt.Legend(title="Cell Type")), # Customize legend title
      tooltip=['Sample', 'Cell_Type', alt.Tooltip('Fraction:Q', format='.4f')] # Tooltip for interactivity, format fraction
  ).properties(
      title='Stacked Area Plot of Cell Type Fractions Across Samples'
  )

  chart.show()



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
