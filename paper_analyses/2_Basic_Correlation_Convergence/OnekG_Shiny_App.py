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
import io

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
        ui.input_slider("effect_size", "Effect Sizes", min=0, max=5, value=(0,5), step=0.1),
        ui.input_slider("case_prob", "Case Proportions", min=0, max=0.5, value=(0,0.5), step=0.05),
        ui.input_selectize("color_var", "Color points by", 
                           choices={"n_causal": "# Causal Variants", "n_traits": "# Traits", "effect_size": "Causal Effect Size", "population": "Population", "trait_type": "Trait Type", "case_prob": "Case Proportion"}, 
                           selected="n_traits"),
        ui.input_selectize("group_cols",
            "Group stats table by:",
            choices=["chromosome", "population", "trait_type", "n_traits", "n_causal", "effect_size", "case_prob"],
            multiple=True,
            selected=["chromosome", "population", "trait_type", "n_traits", "n_causal"]),
        ui.input_action_button("update", "Update Plots and Table", class_="btn-primary", disabled=True),
        width=300
    ),
    ui.page_fluid(
        ui.row(ui.output_plot("corr_plot")),
        ui.row(ui.output_data_frame("stats_table")),
        ui.row(ui.download_button("download_plot", "Download Plot as PNG", class_="btn-secondary"),
               ui.download_button("download_table", "Download Table as CSV", class_="btn-secondary")
               )
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
            ((df["case_prob"].between(*input.case_prob())) | (df['trait_type'] == 'continuous')) &
            (df["maf1"].gt(input.min_maf()) if input.min_maf() else True) &
            (df["maf2"].gt(input.min_maf()) if input.min_maf() else True) &
            (df["effect_size"].between(*input.effect_size()) | (df['n_causal'] == 0))
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
            stats = df.groupby(list(group_cols), observed=False, dropna=False)[["geno_corr", "beta_corr"]].agg(
                MAE=("beta_corr", lambda x: np.nanmean(np.abs(x - df.loc[x.index, "geno_corr"]))),
                RMSE=("beta_corr", lambda x: np.sqrt(np.nanmean((x - df.loc[x.index, "geno_corr"])**2))),
                r=("geno_corr", lambda x: x.corr(df.loc[x.index, "beta_corr"]))
            ).reset_index()
            stats = stats[stats['MAE'].notna() & stats['RMSE'].notna() & stats['r'].notna()]
        return stats.round(4)
    
    @render.download(filename="correlation_plot.png")
    def download_plot():
        df = filtered_data()
        if df.empty:
            # Return an empty image if no data
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No data to plot", ha='center', va='center')
        else:
            plot_df = df.sample(n=min(10000, len(df)), random_state=42) if len(df) > 10000 else df
            color_var = input.color_var()
            fig, ax = plt.subplots(figsize=(10, 8))
            # (Repeat your plotting code here, or refactor plotting into a helper function)
            if plot_df[color_var].dtype.kind in "O" or plot_df[color_var].dtype.name == "category":
                categories = plot_df[color_var].astype(str).unique()
                colors = plt.cm.tab10.colors
                color_map = {cat: colors[i % len(colors)] for i, cat in enumerate(categories)}
                point_colors = plot_df[color_var].astype(str).map(color_map)
                scatter = ax.scatter(plot_df["geno_corr"], plot_df["beta_corr"], c=point_colors, alpha=0.5, s=20)
                for cat in categories:
                    ax.scatter([], [], color=color_map[cat], label=cat)
                ax.legend(title=color_var)
            else:
                scatter = ax.scatter(
                    plot_df["geno_corr"],
                    plot_df["beta_corr"],
                    c=plot_df[color_var],
                    cmap="viridis",
                    alpha=0.5,
                    s=20,
                )
                plt.colorbar(scatter, ax=ax, label=color_var)
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
        with io.BytesIO() as buf:
            fig.savefig(buf, format="png")
            plt.close(fig)
            yield buf.getvalue()
    
    @render.download(filename="stats_table.csv")
    def download_table():
        df = filtered_data()
        group_cols = input.group_cols()
        if not group_cols:
            stats = pd.DataFrame({
                "MAE": [np.mean(np.abs(df["beta_corr"] - df["geno_corr"]))],
                "RMSE": [np.sqrt(np.mean((df["beta_corr"] - df["geno_corr"])**2))],
                "r": [df["geno_corr"].corr(df["beta_corr"])]
            })
        else:
            stats = df.groupby(list(group_cols), observed=False, dropna=False)[["geno_corr", "beta_corr"]].agg(
                MAE=("beta_corr", lambda x: np.nanmean(np.abs(x - df.loc[x.index, "geno_corr"]))),
                RMSE=("beta_corr", lambda x: np.sqrt(np.nanmean((x - df.loc[x.index, "geno_corr"])**2))),
                r=("geno_corr", lambda x: x.corr(df.loc[x.index, "beta_corr"]))
            ).reset_index()
            stats = stats[stats['MAE'].notna() & stats['RMSE'].notna() & stats['r'].notna()]
        with io.StringIO() as buf:
            stats.round(4).to_csv(buf, index=False)
            yield buf.getvalue()

app = App(app_ui, server)
