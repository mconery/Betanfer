'''
OnekG_Shiny_App.py

The purpose of this script is to provide a way to easily view the results of 
the 1000G simulations and compare different parameter settings.

***Run in Conda finemap environment***
'''

from shiny import App, render, ui, reactive
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

# Load combined data (modify path to your CSV file)
DATA_PATH = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Data/1000G_Simulations/master.results.csv.gz"

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize("chrom", "Chromosome", choices=[], multiple=True),
        ui.input_selectize("pop", "Population", choices=[], multiple=True),
        ui.input_selectize("trait_type", "Trait Type", choices=[], multiple=True),
        ui.input_selectize("n_traits", "Number of Simulated Traits", choices=[], multiple=True),
        ui.input_slider("min_maf", "Minimum Minor Allele Frequency", min=0, max=0.5, value=0),
        ui.input_slider("n_causal", "Causal Variants", min=0, max=5, value=(0,5)),
        ui.input_selectize("effect_size", "Effect Sizes", choices=[], multiple=True),
        ui.input_selectize("color_var", "Color points by", 
                           choices={"n_causal": "# Causal Variants", "n_traits": "# Traits", "effect_size": "Causal Effect Size", "population": "Population", "trait_type": "Trait Type"}, 
                           selected="n_traits"),
        ui.input_selectize("group_cols",
            "Group stats table by:",
            choices=["chromosome", "population", "trait_type", "n_traits", "n_causal", "effect_size"],
            multiple=True,
            selected=["chromosome", "population", "trait_type", "n_traits", "n_causal"]),
        ui.input_action_button("update", "Update Plots and Table", class_="btn-primary", disabled=True),
        width=300
    ),
    ui.page_fluid(
        ui.row(ui.output_plot("corr_plot")),
        ui.row(ui.output_data_frame("stats_table"))
    )
)

def server(input, output, session):
    @reactive.calc
    def loaded_data():
        temp = pd.read_csv(
            DATA_PATH,
            sep=',',
            dtype={'chromosome': str, 'population': str, 'trait_type': str, 'n_traits': str}
        )
        temp['n_traits'] = pd.Categorical(temp['n_traits'], categories=['50', '100', '200', '300', '400', '500', '1000'])
        return temp
    
    @reactive.effect
    def update_inputs():
        df = loaded_data()
        ui.update_selectize("chrom", choices=df["chromosome"].astype(str).unique().tolist())
        ui.update_selectize("pop", choices=df["population"].astype(str).unique().tolist())
        ui.update_selectize("trait_type", choices=df["trait_type"].astype(str).unique().tolist())
        ui.update_selectize("n_traits", choices=df["n_traits"].astype(str).unique().tolist()) 
        mafs = np.array([0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
        ui.update_selectize("min_maf", choices=[str(x) for x in mafs])
        effect_sizes = df["effect_size"].dropna().round(2).unique()
        ui.update_selectize("effect_size", choices=[str(x) for x in effect_sizes])
        ui.update_action_button("update", disabled=False)

    @reactive.calc
    def filtered_data():
        df = loaded_data()
        return df[
            (df["chromosome"].isin(input.chrom()) if input.chrom() else True) &
            (df["population"].isin(input.pop()) if input.pop() else True) &
            (df["trait_type"].isin(input.trait_type()) if input.trait_type() else True) &
            (df["n_traits"].astype(str).isin(input.n_traits()) if input.n_traits() else True) & 
            df["n_causal"].between(*input.n_causal()) &
            (df["maf1"].gt(input.min_maf()) if input.min_maf() else True) &
            (df["maf2"].gt(input.min_maf()) if input.min_maf() else True) &
            (df["effect_size"].isin(input.effect_size()) if input.effect_size() else True)
        ]

    @render.plot(alt="Beta vs Genotype Correlations")
    @reactive.event(input.update)
    def corr_plot():
        df = filtered_data()
        if df.empty:
            return plt.figure()
    
        plot_df = df.sample(n=min(10000, len(df)), random_state=42) if len(df) > 10000 else df
    
        color_var = input.color_var()
        fig, ax = plt.subplots(figsize=(10, 8))
    
        # Handle categorical and numeric coloring
        if plot_df[color_var].dtype.kind in "O" or plot_df[color_var].dtype.name == "category":
            # For categorical variables, assign a color to each category
            categories = plot_df[color_var].astype(str).unique()
            colors = plt.cm.tab10.colors  # Or another categorical colormap
            color_map = {cat: colors[i % len(colors)] for i, cat in enumerate(categories)}
            point_colors = plot_df[color_var].astype(str).map(color_map)
            scatter = ax.scatter(
                plot_df["geno_corr"], plot_df["beta_corr"], c=point_colors, alpha=0.5, s=20
            )
            # Create legend for categories
            for cat in categories:
                ax.scatter([], [], color=color_map[cat], label=cat)
            ax.legend(title=color_var)
        else:
            # For numeric variables, use a continuous colormap
            scatter = ax.scatter(
                plot_df["geno_corr"],
                plot_df["beta_corr"],
                c=plot_df[color_var],
                cmap="viridis",
                alpha=0.5,
                s=20
            )
            plt.colorbar(scatter, ax=ax, label=color_var)
    
        # Add LOWESS trend line if desired
        try:
            if len(plot_df) > 100:
                lowess_line = lowess(plot_df["beta_corr"], plot_df["geno_corr"], frac=0.2)
                ax.plot(lowess_line[:, 0], lowess_line[:, 1], color='red', linewidth=2)
        except ImportError:
            pass
    
        ax.set_xlabel("Genotype Correlation")
        ax.set_ylabel("Beta Correlation")
        ax.set_title("Beta vs Genotype Correlations")
        plt.tight_layout()
        return fig


    @render.data_frame
    @reactive.event(input.update)
    def stats_table():
        df = filtered_data()
        group_cols = input.group_cols()
        if not group_cols:
            # If no grouping columns selected, summarize the entire filtered dataset
            stats = pd.DataFrame({
                "MAE": [np.mean(np.abs(df["beta_corr"] - df["geno_corr"]))],
                "RMSE": [np.sqrt(np.mean((df["beta_corr"] - df["geno_corr"])**2))],
                "r": [df["geno_corr"].corr(df["beta_corr"])]
            })
        else:
            stats = df.groupby(list(group_cols), observed=False)[["geno_corr", "beta_corr"]].agg(
                MAE=("beta_corr", lambda x: np.mean(np.abs(x - df.loc[x.index, "geno_corr"]))),
                RMSE=("beta_corr", lambda x: np.sqrt(np.mean((x - df.loc[x.index, "geno_corr"])**2))),
                r=("geno_corr", lambda x: x.corr(df.loc[x.index, "beta_corr"]))
            ).reset_index()
        return stats.round(4)


app = App(app_ui, server)
