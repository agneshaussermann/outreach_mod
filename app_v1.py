from shiny import App, ui, render, reactive
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import io
import base64
import sys


## Build model

import numpy
import cobra
import pandas
from cobra import Model, Reaction, Metabolite

# metabolites

co2 = Metabolite(
    'co2',
    name='Carbon dioxide',
    compartment='e0')

gluc = Metabolite(
    'glc__D',
    name='D-glucose',
    compartment='e0')

n2 = Metabolite(
    'n2',
    name='Nitrogen',
    compartment='e0')

nh4 = Metabolite(
    'nh4',
    name='Nutrients',
    compartment='e0')

o2 = Metabolite(
    'o2',
    name='O2',
    compartment='e0')

bio = Metabolite(
    'biomass',
    name='Total biomass',
    compartment='e0')

aut_bio = Metabolite(
    'auto_biomass',
    name='Autotrophic biomass',
    compartment='e0')

het_bio = Metabolite(
    'het_biomass',
    name='Hetrotrophic biomass',
    compartment='e0')

diaz_bio = Metabolite(
    'diaz_biomass',
    name='Diazotrophic biomass',
    compartment='e0')

hn = Metabolite(
    'photon',
    name = 'Photon',
    compartment = 'e0')

# model - build new model for every input

def build_model(auto, hetero, diaz, ecosystem):
    
    tot_abun = auto + hetero + diaz
    auto_rel = auto / tot_abun
    het_rel = hetero / tot_abun
    diaz_rel = diaz / tot_abun

    model = Model("community")
    
    reaction1 = Reaction('autotroph')
    reaction1.name = 'Carbon dioxide fixation'
    reaction1.lower_bound = 0. 
    reaction1.upper_bound = auto  # scaled to abundance

    reaction1.add_metabolites({
        co2: -1.0,
        hn: -1.0,
        nh4: -1.0,
        gluc: 1.0,
        o2: 1.0,
        aut_bio: 1

    })

    reaction2 = Reaction('heterotroph')
    reaction2.name = 'Glucose oxidation'
    reaction2.lower_bound = 0.  
    reaction2.upper_bound = hetero 

    reaction2.add_metabolites({
        co2: 4.0,
        gluc: - 2.0,
        nh4: - 1.0,
        o2: -1.0,
        het_bio: 1.0
    })

    reaction3 = Reaction('diazotroph')
    reaction3.name = 'Nitrogen fixation'
    reaction3.lower_bound = 0.  
    reaction3.upper_bound = diaz 

    reaction3.add_metabolites({
        nh4: 1.0,
        co2: 1.0,
        n2: - 1.0,
        gluc: -1.0,
        diaz_bio: 1.0
    })


    biomass_rxn = Reaction('biomass')
    biomass_rxn.name = 'Total biomass formation'
    biomass_rxn.lower_bound = 0.  # This is the default
    biomass_rxn.upper_bound = 1000.  # This is the default

    biomass_rxn.add_metabolites({
        aut_bio: - auto_rel,
        het_bio: - het_rel,
        diaz_bio: - diaz_rel,
        bio: 1.0
    })

#add metabolites
    model.add_metabolites([co2, gluc, bio, hn, aut_bio, het_bio, diaz_bio, nh4, n2, o2])

#add reactions
    model.add_reactions([reaction1, reaction2, reaction3, biomass_rxn])

#add exchange reactions
    model.add_boundary(model.metabolites.get_by_id("co2"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("photon"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("glc__D"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("nh4"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("n2"), type="exchange")
    #model.add_boundary(model.metabolites.get_by_id("no3"), type="exchange")
    #model.add_boundary(model.metabolites.get_by_id("co2"), type="sink")
    #model.add_boundary(model.metabolites.get_by_id("glc__D"), type="sink")
    model.add_boundary(model.metabolites.get_by_id("o2"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("biomass"), type="exchange")
    
#medium ie ecosystem conditions

    if ecosystem == "lake":
        medium = model.medium
        medium["EX_co2"] = 10.0 # change for open system
        medium["EX_photon"] = 100 #light
        medium["EX_glc__D"] = 1.0 # change for open system
        medium["EX_n2"] = 100
        medium["EX_o2"] = 1.0
        medium["EX_nh4"] = 1.0
        medium["EX_biomass"] = 0 
        model.medium = medium
    
    model.objective = "biomass"
    
    solution = model.optimize()
    
    return solution


## animation + shiny ui

from matplotlib.patches import FancyArrowPatch

# =========================================================
# DRAWING UTILITIES
# =========================================================
def draw_arrow(ax, p1, p2, width, color, rad=0.15):
    ax.add_patch(FancyArrowPatch(
        p1, p2,
        arrowstyle="-|>",
        connectionstyle=f"arc3,rad={rad}",
        linewidth=width,
        color=color,
        alpha=0.85,
        mutation_scale=18
    ))

import plotly.graph_objects as go
from shinywidgets import output_widget


def interp(p1, p2, t):
    return (
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
    )


def add_particle_stream(
    fig,
    p1,
    p2,
    color,
    phase,
    n_particles=4,
    size=10
):
    xs = []
    ys = []

    for k in range(n_particles):

        t = (phase + k / n_particles) % 1

        x, y = interp(p1, p2, t)

        xs.append(x)
        ys.append(y)

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers",
            marker=dict(
                size=size,
                color=color
            ),
            showlegend=False,
            hoverinfo="skip"
        )
    )


def draw_frame(frame, sol):

    flux = sol.fluxes
    flux_abs = flux.abs()

    pos = {
        "Light": (0.1, 0.9),
        "CO2": (0.45, 0.75),
        "N2": (0.8, 0.9),

        "Autotrophs": (0.25, 0.65),
        "Heterotrophs": (0.15, 0.35),

        "Diazotrophs": (0.65, 0.65),
        "NH4": (0.6, 0.45),

        "Carbohydrates": (0.35, 0.25),
    }

    fig = go.Figure()

    # --------------------------------------------------
    # Static network
    # --------------------------------------------------

    edges = [
        ("Light", "Autotrophs", "green", "autotroph"),
        ("CO2", "Autotrophs", "green", "autotroph"),
        ("N2", "Diazotrophs", "purple", "diazotroph"),

        ("Diazotrophs", "NH4", "purple", "diazotroph"),

        ("NH4", "Autotrophs", "green", "autotroph"),
        ("NH4", "Heterotrophs", "red", "heterotroph"),

        ("Autotrophs", "Carbohydrates", "green", "autotroph"),
        ("Carbohydrates", "Heterotrophs", "red", "heterotroph"),

        ("Heterotrophs", "CO2", "blue", "heterotroph"),
    ]

    max_flux = max(flux_abs.max(), 1e-6)

    def scale_flux(f):
        #return 1 + 10 * (f / max_flux)
        return 10 * f

    for src, dst, color, rxn in edges:

        x0, y0 = pos[src]
        x1, y1 = pos[dst]

        width = scale_flux(flux_abs[rxn])

        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(
                    color=color,
                    width=width
                ),
                opacity=0.5,
                showlegend=False,
                hoverinfo="skip"
            )
        )

    # --------------------------------------------------
    # Nodes
    # --------------------------------------------------

    fig.add_trace(
        go.Scatter(
            x=[p[0] for p in pos.values()],
            y=[p[1] for p in pos.values()],
            text=list(pos.keys()),
            mode="markers+text",
            marker=dict(
                size=45,
                color="lightblue",
                line=dict(color="black", width=1)
            ),
            textposition="middle center",
            showlegend=False
        )
    )

    # --------------------------------------------------
    # Timeline
    # --------------------------------------------------

    if frame < 20:

        phase = frame / 20

        add_particle_stream(
            fig,
            pos["Light"],
            pos["Autotrophs"],
            "green",
            phase,
        )

        add_particle_stream(
            fig,
            pos["N2"],
            pos["Diazotrophs"],
            "purple",
            phase,
        )

    elif frame < 50:

        phase = (frame - 20) / 30

        add_particle_stream(
            fig,
            pos["Autotrophs"],
            pos["Carbohydrates"],
            "green",
            phase,
        )

        add_particle_stream(
            fig,
            pos["Diazotrophs"],
            pos["NH4"],
            "purple",
            phase,
        )

        add_particle_stream(
            fig,
            pos["NH4"],
            pos["Autotrophs"],
            "green",
            phase,
        )

    elif frame < 90:

        phase = (frame - 70) / 20

        add_particle_stream(
            fig,
            pos["Carbohydrates"],
            pos["Heterotrophs"],
            "red",
            phase,
        )

        add_particle_stream(
            fig,
            pos["NH4"],
            pos["Heterotrophs"],
            "red",
            phase,
        )

    else:

        phase = (frame - 90) / 10

        add_particle_stream(
            fig,
            pos["Heterotrophs"],
            pos["CO2"],
            "blue",
            phase,
        )

    fig.update_layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white",
        margin=dict(l=10, r=10, t=10, b=10),
        height=600,
    )

    return fig

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_numeric("auto", "Autotroph flux", 1),
        ui.input_numeric("hetero", "Heterotroph flux", 1),
        ui.input_numeric("diaz", "Diazotroph flux", 1),
        ui.input_select("ecosystem", "Ecosystem", ["lake"]),
        ui.input_slider("frame", "Frame", 0, 100, 0),
    ),
    output_widget("anim"),
    ui.tags.script("""
        let f = 0;
        setInterval(() => {
            f = (f + 1) % 100;
            Shiny.input.set("frame", f, {priority: "event"});
        }, 100);
    """),
    title="Flux Animation"
)

# ---------------------------------------------------------
# Server
# ---------------------------------------------------------
from shinywidgets import render_plotly

def server(input, output, session):

    @reactive.calc
    def solution():
        return build_model(
            input.auto(),
            input.hetero(),
            input.diaz(),
            input.ecosystem()
        )

    @output
    @render_plotly
    def anim():
        return draw_frame(
            input.frame(),
            solution()
        )


app = App(app_ui, server)
