#!/bin/bash
# EpiTox revision — GitHub issues
# Usage: gh auth login (if not already), then bash create_epitox_issues.sh
# Replace REPO with your repo e.g. hoor-git/EpiTox or org/repo

REPO="HA1-biocopy/EpiTox_paper"

# ── Create labels first ────────────────────────────────────────────────────
gh label create "priority:high"   --color "B60205" --description "Must address before submission"       --repo "$REPO" 2>/dev/null
gh label create "priority:medium" --color "E4A11B" --description "Important but not blocking"           --repo "$REPO" 2>/dev/null
gh label create "priority:low"    --color "0E8A16" --description "Minor fix or polish"                  --repo "$REPO" 2>/dev/null
gh label create "pending"         --color "D4C5F9" --description "Waiting on external input"            --repo "$REPO" 2>/dev/null
gh label create "frontiers"       --color "1D76DB" --description "Specifically important for Frontiers" --repo "$REPO" 2>/dev/null

# ── NEW SCIENCE ────────────────────────────────────────────────────────────

gh issue create \
  --repo "$REPO" \
  --title "Write SNP-based answer to NetMHCpan" \
  --label "priority:high,frontiers" \
  --body "## What
Write the core argument that directly answers the reviewers' NetMHCpan objection.

## Why
Reviewers argue: affinity is the deciding factor, NetMHCpan gives affinity, so why use EpiTox?

The answer is categorical, not a performance comparison:
- EpiTox is not an HLA-binding predictor — NetMHCpan is a component inside EpiTox, not a competitor
- SNP-derived peptides do not exist in the reference proteome — NetMHCpan cannot generate or score them
- 2 of 3 confirmed SNP binders exceeded target peptide affinity — invisible to any affinity-only approach by design

## Location
Main paper — Discussion section

## Effort
~3–4 days"

gh issue create \
  --repo "$REPO" \
  --title "Recall/precision benchmark: FindTopes vs Expitope and iCrossR" \
  --label "priority:high,frontiers" \
  --body "## What
Add a quantitative benchmark table comparing FindTopes against existing tools on a common ground truth dataset.

## Why
Frontiers reviewers expect benchmark numbers, not just feature comparison tables. Reviewer 1 and the Associate Editor both flagged the absence of quantitative comparison explicitly.

## Approach
- Use the 19 known MAGEA3 off-targets as ground truth
- Compute recall, precision, and F1 for FindTopes, Expitope, and iCrossR on the same set
- Decide location (main vs supplementary) based on space and narrative flow

## Effort
~1 week"

gh issue create \
  --repo "$REPO" \
  --title "Quantitative benchmark vs affinity-only baseline" \
  --label "priority:high" \
  --body "## What
Show that EpiTox's prioritization adds value beyond simply applying a NetMHCpan affinity threshold.

## Why
The 9.7% hit rate is difficult to interpret without a baseline. Reviewers need to see enrichment, not a raw number.

## Approach
- Define the affinity-passing pool (all peptides with KD ≤ 500 nM)
- Compare hit rate across EpiTox risk tiers vs random draw from that pool
- Acknowledge explicitly that affinity alone retrieves the 36 binders because affinity was a selection criterion — circular by design, not a failure

## Location
Main paper or supplementary — TBD based on space

## Effort
~1 week"

gh issue create \
  --repo "$REPO" \
  --title "[PENDING] Cell assay data — full 538 peptide set" \
  --label "priority:high,pending" \
  --body "## What
Incorporate cell loading / functional assay data for the full 538 peptide validation set if available.

## Why
Cell context data would directly answer Reviewer 2's major criticism that binding ≠ downstream biological effects. This is potentially one of the strongest additions to the paper.

## Status
Cell assays are being established. Data expected within ~1 week. **Do not start writing until data is confirmed available.**

## Location
Main paper — Results and Discussion

## Effort
TBD depending on data completeness"

# ── SCIENTIFIC CLARIFICATIONS ─────────────────────────────────────────────

gh issue create \
  --repo "$REPO" \
  --title "Selection bias clarification — bi-feature score validation history" \
  --label "priority:high" \
  --body "## What
Explicitly clarify the validation pool selection strategy and its relationship to the scoring timeline.

## Why
Reviewers flagged that the tested pool skews heavily toward high bi-feature risk peptides (~80%), suggesting inflated hit rates. The honest explanation is stronger than an apology.

## What to write
At the time of experimental validation, only the bi-feature score existed. The pool was selected using that score — this is why the pool skews high-risk. The multi-feature score and evidence-based framework were developed subsequently. The paper demonstrates their complementarity retrospectively, which is a methodological strength, not a bias.

## Location
Methods section (validation candidate selection) + Results framing

## Effort
~1 day"

gh issue create \
  --repo "$REPO" \
  --title "Justify multi-feature score design — why these features, why target-aware strategy" \
  --label "priority:medium" \
  --body "## What
Add explicit scientific justification for the multi-feature score design choices.

## Why
Reviewer 2 flagged that the physicochemical feature selection is 'guided by observations from internal analyses' without external justification. The target-aware weighting strategy also needs motivation.

## What to address
- Why these 8 physicochemical features specifically?
- Why backbone vs global feature distinction based on anchor status?
- Why is target-aware strategy (weighting by mismatch location relative to anchors) the right approach here vs a generic similarity metric?
- If internal datasets informed the choice, can any of that rationale be made explicit or supported by external references?

## Location
Methods section — AnnoTopes / multi-feature score

## Effort
~1–2 days"

gh issue create \
  --repo "$REPO" \
  --title "Explain single-antigen validation depth as deliberate methodological choice" \
  --label "priority:medium" \
  --body "## What
State explicitly why MAGEA3 was used as the sole validation target.

## Why
All three reviewers flagged the single antigen as a limitation. The counter-argument needs to be written directly, not left implicit.

## Argument to make
MAGEA3 was selected because it provides 19 experimentally validated off-targets as ground truth — a requirement for rigorous benchmarking. Validation depth (538 peptides, quantitative KD measurements) over breadth is an intentional methodological choice. Other tools use multiple antigens precisely because they cannot validate at this scale. MAGEA4 (ref 28) provides independent corroboration of generalizability.

## Location
Methods rationale paragraph + Discussion

## Effort
~half day"

# ── FRAMING & WRITING ─────────────────────────────────────────────────────

gh issue create \
  --repo "$REPO" \
  --title "Frame EpiTox as methods + validation paper, not a software paper" \
  --label "priority:high" \
  --body "## What
Adjust framing throughout to position EpiTox as a methods and experimental validation paper rather than a software distribution paper.

## Why
Most bioinformatics journals' mandatory code-sharing policies specifically target software submissions. EpiTox's primary contribution is the computational framework combined with 538-peptide experimental validation. This framing distinction matters for the repo question and should be reflected in the abstract, introduction, and cover letter.

## Where to apply
- Abstract: lead with framework + experimental validation, not tool description
- Introduction: frame contribution as methodology + validation, not software release
- Cover letter: explicitly state this distinction

## Effort
~1 day"

gh issue create \
  --repo "$REPO" \
  --title "Highlight MAGEA4 paper (ref 28) in multiple places as corroborating evidence" \
  --label "priority:medium" \
  --body "## What
Reference the published MAGEA4 paper in multiple strategic locations throughout the manuscript.

## Why
Currently ref 28 appears as a single footnote in the Discussion. The meeting agreed it should appear in multiple places as corroborating evidence for generalizability — but not as a second case study, which would confuse the MAGEA3 storyline.

## Where to add references
- Introduction: when establishing the broader applicability of the approach
- Results: where generalizability is implied
- Discussion: explicitly state that independent validation on MAGEA4 confirms consistent off-target identification

## Effort
~1 day"

gh issue create \
  --repo "$REPO" \
  --title "Strengthen practical utility framing for non-specialist readers" \
  --label "priority:medium,frontiers" \
  --body "## What
Sharpen the framing of EpiTox as a practical decision-support tool accessible to non-computational drug safety scientists.

## Why
Frontiers in Bioinformatics explicitly values guidance to non-specialists. The decision matrix and modular architecture are genuine strengths here but are currently underemphasized.

## What to improve
- Introduction: frame EpiTox as a toolkit a drug safety scientist can follow, not just a bioinformatics tool
- Decision matrix (Table 1): make the experimental goal → approach mapping more prominent
- Methods: ensure module descriptions are accessible without deep computational background

## Effort
~1 day"

# ── MINOR FIXES ───────────────────────────────────────────────────────────

gh issue create \
  --repo "$REPO" \
  --title "Minor text fixes — definitions, LR justification, figure clarity" \
  --label "priority:low" \
  --body "## What
Address the minor but specific points raised by Reviewer 2.

## Checklist
- [ ] Define PEPREP before first use in main text (currently appears without introduction)
- [ ] Define 'backbone' and 'exposed positions' explicitly before they are used for classification
- [ ] Justify LR values with external references where possible, not just internal sensitivity analysis
- [ ] Clarify the 538 peptide selection criteria more explicitly in Methods
- [ ] Improve Fig 1 readability — consider simplifying or splitting into panels

## Effort
~2 days"

echo ""
echo "Done. $(gh issue list --repo "$REPO" --limit 20 | wc -l) issues created in $REPO"
