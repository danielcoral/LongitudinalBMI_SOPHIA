"""Analysis dependencies; publication exports are explicit opt-in targets."""
from pathlib import Path
import json
import os
import re

configfile: "config/config.yaml"

stages = config["stages"]
paths = config["paths"]
helpers = sorted(str(path) for path in Path("src").glob("*.R"))
main_figures = json.loads(Path("paper/mainfigsmap.json").read_text())
supp_tables = json.loads(Path("paper/supptabsmap.json").read_text())
supp_images = [
  os.path.normpath(str(Path("paper") / image))
  for image in re.findall(
    r"!\[[^\]]*\]\(([^)]+)\)", Path("paper/suppfigsmap.md").read_text()
  )
]


def resolved(files):
  return [
    path.format(**paths, **config["cohort_results"]) for path in files
  ]


rule all:
  input:
    sorted({
      path for stage in stages.values()
      for path in resolved(stage["outputs"])
    })


for stage_name, stage in stages.items():
  stage_data = resolved(stage["inputs"])
  if stage_name == "stage_29":
    stage_data.append("paper/supp_notes.md")
  rule:
    name: stage_name
    input:
      script=stage["script"],
      data=stage_data,
      code=helpers + ["config/config.yaml", "environment.yml"]
    output:
      resolved(stage["outputs"])
    threads:
      stage["threads"]
    resources:
      mem_mb=stage["mem_mb"]
    log:
      f"logfiles/{stage_name}.log"
    shell:
      "Rscript --vanilla {input.script:q} > {log:q} 2>&1"


rule diagrams:
  input:
    f"{paths['plots']}/study_flowchart.png",
    f"{paths['plots']}/graphical_abstract.png"


for number, stem in [("25", "study_flowchart"), ("26", "graphical_abstract")]:
  rule:
    name: f"diagram_{number}"
    input:
      script=next(str(path) for path in Path("scripts").glob(f"{number}_*.py")),
      settings=["environment.yml", "config/config.yaml"]
    output:
      f"{paths['plots']}/{stem}.png"
    log:
      f"logfiles/diagram_{number}.log"
    shell:
      "python {input.script:q} > {log:q} 2>&1"


rule manuscript:
  input:
    "paper/manuscript.md", "paper/references.bib",
    "paper/formats/vancouver-brackets.csl", "paper/formats/reference.docx",
    "paper/scripts/buildmanus.sh", "environment.yml"
  output:
    "paper/manuscript.docx"
  shell:
    "bash paper/scripts/buildmanus.sh"


rule supplementary_figures:
  input:
    "paper/suppfigsmap.md", "paper/formats/reference.docx",
    "paper/scripts/pagebreak.lua", "paper/scripts/buildsuppfigs.sh",
    "paper/scripts/formatsuppfigs.py", "environment.yml", supp_images
  output:
    "paper/supp_figs.docx"
  shell:
    "bash paper/scripts/buildsuppfigs.sh"


rule supplementary_notes:
  input:
    "paper/supp_notes.md", "paper/references.bib",
    "paper/formats/vancouver-brackets.csl", "paper/formats/reference.docx",
    "paper/scripts/buildsuppnotes.sh", "environment.yml"
  output:
    "paper/supp_notes.docx"
  shell:
    "bash paper/scripts/buildsuppnotes.sh"


rule outcome_figures:
  input:
    "paper/scripts/buildoutcomesfigs.R", "src/project_setup.R",
    "config/config.yaml", "environment.yml",
    resolved([
      "{ukb_results}/table_incident_outcomes_marginal_quintiles.csv",
      "{ukb_processed}/plot_adjusted_outcomes.rds",
      "{ukb_processed}/plot_dpp_comparison.rds"
    ])
  output:
    resolved(["{plots}/incident_outcomes_summary.png",
              "{plots}/dpp_comparison.png"])
  shell:
    "Rscript --vanilla paper/scripts/buildoutcomesfigs.R"


rule main_figures:
  input:
    "paper/mainfigsmap.json", "paper/scripts/buildmainfigs.R",
    "environment.yml",
    [entry["file"] for entry in main_figures.values()],
    [entry["source"] for entry in main_figures.values()]
  output:
    [
      f"paper/mainfigs/{name}{Path(entry['file']).suffix}"
      for name, entry in main_figures.items()
    ]
  shell:
    "Rscript --vanilla paper/scripts/buildmainfigs.R"


rule supplementary_tables:
  input:
    "paper/supptabsmap.json", "paper/scripts/buildsupptabs.R",
    "environment.yml",
    [entry["file"] for entry in supp_tables.values()],
    [entry["source"] for entry in supp_tables.values()]
  output:
    "paper/supp_tables.xlsx"
  shell:
    "Rscript --vanilla paper/scripts/buildsupptabs.R"
