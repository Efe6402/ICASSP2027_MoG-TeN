NETWORK REDUCTION — THREE-EXPERIMENT NUCLEAR-ABLATION SUITE
================================================================

This package contains the final three synthetic experiments for the paper.

The SAME nuclear-norm ablation is now included as a baseline in all three:
    Proposed
    Proposed (tau=0)
    BTD
    CPD

The tau=0 baseline is PAIRED:
after the positive-tau Proposed hyperparameters are selected, ONLY tau is
changed to zero. gamma, mu, lambda, and zeta remain fixed.

COMMON HYPERPARAMETER SELECTION
-------------------------------
OURS hyperparameters are selected using CLEAN component recovery, not clean-X
reconstruction:

    maximize min(mean clean-C correlation,
                 mean clean-alpha correlation).

Clean-X NMSE is used only for evaluation.

RUNME 1
-------
RUNME_1_noise_sweep_with_tau_ablation.m

q=2, eta sweep:
    0.316, 0.422, 0.562, 0.750, 1, 1.334, 1.778, 2.239

Submission output:
    submission_noise_sweep_X_NMSE.png

Diagnostics:
    C_r correlation vs eta
    alpha_r correlation vs eta
    balanced factor score vs eta


RUNME 2
-------
RUNME_2_q_sweep_with_tau_ablation.m

eta = 0.316227766 fixed
q = 2,4,8,16

Submission output:
    submission_q_sweep_X_NMSE.png

Diagnostics:
    C_r correlation vs q
    alpha_r correlation vs q
    balanced factor score vs q


RUNME 3
-------
RUNME_3_sigma_sweep_with_tau_ablation.m

RESWEEP ONLY
q = 2
eta = 0.5623413252
sigma = 1.5,2,2.5,3,3.5,4,4.5,5

Submission output:
    submission_sigma_sweep_X_NMSE.png

Diagnostics:
    C_r correlation vs sigma
    alpha_r correlation vs sigma
    balanced factor score vs sigma


SUBMISSION FIGURES
------------------
Each RUNME saves its primary figure inside its timestamped results folder AND
copies the latest version into:

    submission_figures/

After running all three scripts, that folder contains the three paper figures.

Primary figures:
    - have NO title (subfigure captions should carry the title)
    - use a simple Tensor X Reconstruction NMSE y-axis
    - use clear academic legends:
          Proposed
          Proposed (tau=0)
          BTD
          CPD
    - use sparse logarithmic eta ticks for the noise sweep
    - save .png (300 dpi), .pdf vector, and .fig


STYLE TUNING
------------
Global submission styling is centralized in:

    lib/submission_axes.m
    lib/submission_labels.m
    lib/submission_legend.m

Comments in those files show exactly where to change:
    font family
    tick font size
    axis-label font size
    legend font size
    axis thickness

Tensorlab path
--------------
The three RUNMEs use tensorlab 

Change CFG.tensorlabPath if necessary.
