# Publications

Use this page to find papers that are directly connected to SpaceAGORA.jl.
The project citation metadata remains in the repository `CITATION.cff` file;
the papers below document SpaceAGORA-based analyses, verification, and guidance
studies.

## Direct SpaceAGORA papers

- Evan Yu, Nihal Simha, and Giusy Falcone, "A Comparative Analysis of
  Aerobraking Dynamics Across Celestial Bodies Using SpaceAGORA," AIAA SCITECH
  2026 Forum, Orlando, FL, 2026. DOI: [10.2514/6.2026-1692](https://doi.org/10.2514/6.2026-1692).

  This paper uses SpaceAGORA.jl for comparative aerobraking perturbation
  analysis across Mars, Venus, Earth, and Titan, including atmosphere, gravity,
  third-body, and solar-radiation-pressure effects.

- Evan Yu, James Cutler, and Giusy Falcone, "Verification of SpaceAGORA.jl
  6-DoF Dynamics Using University of Michigan SmallSats Telemetry," 2026 IEEE
  Aerospace Conference, Big Sky, MT, 2026, pp. 1-14. DOI:
  [10.1109/AERO66936.2026.11520056](https://doi.org/10.1109/AERO66936.2026.11520056).

  This paper validates SpaceAGORA.jl translational and attitude dynamics
  against University of Michigan SmallSat telemetry and Basilisk comparisons.

## Closely connected implementation papers

- Jakob Gorisek-Gazze and Giusy Falcone, "HyPR: Safety- and Fuel-Aware Motion
  Planning and Replanning for Autonomous SmallSat Servicing," accepted to the
  IEEE International Conference on Space Robotics (iSpaRo), 2026.
  [Lab publication record](https://www.spacefalconlab.com/publications/).

  This study develops a planner for spacecraft motion near obstacles, combining
  an RRT-Connect warm start, particle-swarm refinement of Bézier paths, retiming,
  and LQ-MPC tracking. It reports static and replanning comparisons. The
  [HYPR guide](user/hypr_planning.md) explains the method's connection to
  SpaceAGORA and distinguishes its reported results from the current examples.
  The citation records acceptance; a DOI and proceedings pages are not yet
  established by the cited record.

- Nihal Simha, Evan Yu, and Giusy Falcone, "Generalized Energy Depletion
  Guidance for Aerobraking at Mars, Venus, Earth, and Titan," AIAA SCITECH
  2026 Forum, Orlando, FL, 2026. DOI: [10.2514/6.2026-1378](https://doi.org/10.2514/6.2026-1378).

  This paper implements and exercises the generalized Energy Depletion Guidance
  framework in SpaceAGORA.jl for drag passages at Mars, Venus, Earth, and Titan.
  The [EDG guide](user/edg_guidance.md) describes the current configuration,
  physical quantities and limits of that implementation.

## Foundational guidance paper

- Giusy Falcone and Zachary R. Putnam, "Energy Depletion Guidance for Aerobraking
  Atmospheric Passes," *Journal of Guidance, Control, and Dynamics*, Vol. 45,
  No. 4, 2022, pp. 651–668. DOI:
  [10.2514/1.G006171](https://doi.org/10.2514/1.G006171).

  This paper develops the original EDG formulation for atmospheric passes with
  articulating panels. It provides the foundation for the generalized guidance
  study above; it should not be cited as validation of every later software
  option or adapter.

## BibTeX

```bibtex
@inproceedings{yu2026spaceagoraAerobraking,
  author = {Yu, Evan and Simha, Nihal and Falcone, Giusy},
  title = {A Comparative Analysis of Aerobraking Dynamics Across Celestial Bodies Using SpaceAGORA},
  booktitle = {AIAA SCITECH 2026 Forum},
  address = {Orlando, FL},
  year = {2026},
  doi = {10.2514/6.2026-1692}
}

@inproceedings{yu2026spaceagoraVerification,
  author = {Yu, Evan and Cutler, James and Falcone, Giusy},
  title = {Verification of SpaceAGORA.jl 6-DoF Dynamics Using University of Michigan SmallSats Telemetry},
  booktitle = {2026 IEEE Aerospace Conference},
  address = {Big Sky, MT},
  pages = {1--14},
  year = {2026},
  doi = {10.1109/AERO66936.2026.11520056}
}

@inproceedings{simha2026generalizedEdg,
  author = {Simha, Nihal and Yu, Evan and Falcone, Giusy},
  title = {Generalized Energy Depletion Guidance for Aerobraking at Mars, Venus, Earth, and Titan},
  booktitle = {AIAA SCITECH 2026 Forum},
  address = {Orlando, FL},
  year = {2026},
  doi = {10.2514/6.2026-1378}
}

@inproceedings{gorisekgazze2026hypr,
  author = {Gorisek-Gazze, Jakob and Falcone, Giusy},
  title = {{HyPR}: Safety- and Fuel-Aware Motion Planning and Replanning for Autonomous {SmallSat} Servicing},
  booktitle = {IEEE International Conference on Space Robotics (iSpaRo)},
  year = {2026},
  note = {Accepted},
  url = {https://www.spacefalconlab.com/publications/}
}

@article{falcone2022edg,
  author = {Falcone, Giusy and Putnam, Zachary R.},
  title = {Energy Depletion Guidance for Aerobraking Atmospheric Passes},
  journal = {Journal of Guidance, Control, and Dynamics},
  volume = {45},
  number = {4},
  pages = {651--668},
  year = {2022},
  doi = {10.2514/1.G006171}
}
```
