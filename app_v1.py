from shiny import App, ui, render, reactive
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import io
import base64
import sys
import time


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


# =========================================================
# DRAWING UTILITIES
# =========================================================
def particle_stream(edge_id, color, n_particles):

    out = []

    for k in range(n_particles):

        delay = 4 * k / n_particles

        out.append(f"""
        <circle r="5" fill="{color}">
        </circle>

        <script>
        (() => {{

            const circles =
                document.querySelectorAll(
                    'circle[fill="{color}"]'
                );

            const circle =
                circles[circles.length - 1];

            const path =
                document.getElementById("{edge_id}");

            const length =
                path.getTotalLength();

            function animate(now) {{

                const t =
                    ((now / 1000) + {delay}) % 4;

                const pt =
                    path.getPointAtLength(
                        length * t / 4
                    );

                circle.setAttribute(
                    "cx",
                    pt.x
                );

                circle.setAttribute(
                    "cy",
                    pt.y
                );

                requestAnimationFrame(
                    animate
                );
            }}

            requestAnimationFrame(
                animate
            );

        }})();
        </script>
        """)

    return "\n".join(out)

    
def make_svg(sol):

    flux = sol.fluxes
    flux_abs = flux.abs()

    max_flux = max(flux_abs.max(), 1e-6)

    pos = {
        "Light": (80, 80),
        "CO2": (360, 140),
        "N2": (640, 80),

        "Autotrophs": (220, 220),
        "Diazotrophs": (500, 220),

        "NH4": (450, 360),

        "Carbohydrates": (340, 520),
        "Heterotrophs": (120, 480),
        "O2": (80, 320),
    }

    edges = [
        ("Light", "Autotrophs", "green", "autotroph"),
        ("CO2", "Autotrophs", "green", "autotroph"),

        ("N2", "Diazotrophs", "purple", "diazotroph"),
        ("Diazotrophs", "NH4", "purple", "diazotroph"),

        ("NH4", "Autotrophs", "purple", "autotroph"),
        ("Autotrophs", "Carbohydrates", "green", "autotroph"),

        ("NH4", "Heterotrophs", "purple", "heterotroph"),
        ("Carbohydrates", "Heterotrophs", "red", "heterotroph"),

        ("Heterotrophs", "CO2", "blue", "heterotroph"),

        ("Autotrophs", "O2", "green", "autotroph"),
        ("O2", "Heterotrophs", "red", "heterotroph"),
    ]

    svg = []

    svg.append("""
    <svg id="network" width="800" height="650" viewBox="0 0 800 650">
    """)

    # -----------------------------------
    # Draw paths
    # -----------------------------------

    for i, (src, dst, color, rxn) in enumerate(edges):

        x1, y1 = pos[src]
        x2, y2 = pos[dst]

        width = 12 * flux_abs[rxn]

        svg.append(
            f"""
            <path
                id="edge{i}"
                d="M {x1} {y1} L {x2} {y2}"
                stroke="{color}"
                stroke-width="{width}"
                fill="none"
                opacity="0.6"
            />
            """
        )
    
    # -----------------------------------
    # CO2 exchange arrow
    # -----------------------------------

    ex_flux = flux.get("EX_co2", 0.0)
    print(ex_flux, flux.get("autotroph", 0.0))

    x, y = pos["CO2"]

    if ex_flux > 0:

        w = numpy.exp(ex_flux) + 4
        L = 20 + numpy.exp(ex_flux) * 4

        svg.append(
            f"""
            <polygon
                points="
                {x-w/2},{y-35}
                {x+w/2},{y-35}
                {x+w/2},{y-35-L+15}
                {x+w},{y-35-L+15}
                {x},{y-35-L}
                {x-w},{y-35-L+15}
                {x-w/2},{y-35-L+15}
                "
                fill="darkgreen"
                opacity="0.8"
            />
            """
        )
    # -----------------------------------
    # Nodes
    # -----------------------------------

    for name, (x, y) in pos.items():

        svg.append(
            f"""
            <circle
                cx="{x}"
                cy="{y}"
                r="35"
                fill="lightblue"
                stroke="black"
            />

            <text
                x="{x}"
                y="{y + 5}"
                text-anchor="middle"
                font-size="14"
                font-family="Arial">
                {name}
            </text>
            """
        )

    # -----------------------------------
    # Particles
    # -----------------------------------

    particle_id = 0

    for i, (src, dst, color, rxn) in enumerate(edges):

        n_particles = max(
            0,
            int(
                8 * flux_abs[rxn]
            )
        )

        # sequential pathway activation
        if rxn == "autotroph" and src in ["Light", "CO2"]:
            phase = 0

        elif rxn == "diazotroph":
            phase = 0

        elif src == "NH4":
            phase = 2

        elif rxn == "autotroph" and src in ["Autotrophs"]:
            phase = 2

        elif rxn == "heterotroph" and src != "Heterotrophs":
            phase = 4

        else:
            phase = 6

        for k in range(n_particles):

            pid = f"particle_{particle_id}"
            particle_id += 1

            delay = k / max(n_particles, 1)

            svg.append(
                f"""
                <circle
                    id="{pid}"
                    cx="0"
                    cy="0"
                    r="5"
                    fill="{color}">
                </circle>

                <script>
                (() => {{

                    const particle =
                        document.getElementById("{pid}");

                    const path =
                        document.getElementById("edge{i}");

                    const pathLength =
                        path.getTotalLength();

                    function animate(now) {{

                        const cycle = 10.0;

                        let t =
                            ((now/1000)
                            - {phase}
                            + {delay})
                            % cycle;

                        if (t < 0)
                            t += cycle;

                        const active =
                            (t >= 0 && t <= 2);

                        if (!active) {{

                            particle.setAttribute(
                                "visibility",
                                "hidden"
                            );

                        }} else {{

                            particle.setAttribute(
                                "visibility",
                                "visible"
                            );

                            const frac =
                                t / 2;

                            const point =
                                path.getPointAtLength(
                                    pathLength * frac
                                );

                            particle.setAttribute(
                                "cx",
                                point.x
                            );

                            particle.setAttribute(
                                "cy",
                                point.y
                            );
                        }}

                        requestAnimationFrame(
                            animate
                        );
                    }}

                    requestAnimationFrame(
                        animate
                    );

                }})();
                </script>
                """
            )

    svg.append("</svg>")

    return "\n".join(svg)

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_numeric("auto", "Autotroph flux", 1),
        ui.input_numeric("hetero", "Heterotroph flux", 1),
        ui.input_numeric("diaz", "Diazotroph flux", 1),
        ui.input_select("ecosystem", "Ecosystem", ["lake"]),
    ),
    ui.output_ui("anim"),
    
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
    @render.ui
    def anim():
        return ui.HTML(
            make_svg(
                solution()
            )
        )


app = App(app_ui, server)
