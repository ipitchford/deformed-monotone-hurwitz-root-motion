# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (24-26 November 2026)

1. **Panel of five Claude reviewers** (editor, methodology and proof checking, domain and literature, perspective, and
   devil's advocate) on a merge of three earlier write-ups, with the cited sources read. No mathematical error was
   found and no reviewer raised a critical issue. The major items were editorial: no plain introduction, clashing
   notation between the merged drafts, credit to Accadia, Lewański and Ruzza for the limits, and overclaiming wording.
2. **One revision round with an extension.** A new theorem certifying, for all 132 profiles of degrees 4 to 10,
   real-rootedness at every genus and the exact genera at which zeros move away (Theorem 5.2). Also: a context-first
   introduction and roadmap; one notation; credit to Accadia, Lewański and Ruzza for the limits and for eventual
   real-rootedness and interlacing; Yang's Theorem 1.3 and the Do–He–Robertson conjecture cited; the fixed-point
   identity and Table 1 against CDM's Figure 2; full proofs of the set-partition inversion, the pole multiplicities and
   the occupancy law; scripts converted to explicit checks that stay active under `python -O`.
3. **Re-review** of the revision and **proof audit** of every new or rewritten proof. No mathematical error affecting a
   theorem was found. The items raised (a script that failed under the final replay because of a name clash, an
   attribution between Yang and Do–He–Robertson, a mis-posed open question, four minor proof gaps and the census
   prose against the code) were all corrected before freezing.
4. **Clean replays** of all eight scripts from fresh copies under `python -O`, with a negative control, by the lead and
   again by the publisher (see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)).

No Astra/Codex second-model review took place in this batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| "The limits and eventual real-rootedness already follow from the leading term of Accadia–Lewański–Ruzza." | Agreed and stated in the abstract and Section 1; the new content is the second term and the direction of approach. |
| "The conjecture as stated requires monotonicity from genus 0." | Theorem 5.1 shows it fails from genus 0; the paper proves the eventual form and says so. |
| "The genus threshold is very loose." | Stated after Theorem 4.4; the census shows genus 0 or 2 suffices up to degree 10; Question 5.3 and Section 8. |
| "The fixed-point dependence has an antecedent for simple Hurwitz numbers." | Yang's Theorem 1.3 and the Do–He–Robertson conjecture cited in Section 1 and Remark 3.4. |
| "Several profiles: the direction can be undetermined." | Remark 4.3 gives an explicit case with B = 0 and lists it as open. |
| "Proofs of the inversion, pole multiplicities and occupancy law were compressed." | Written out as Lemmas 2.2-2.3 and the proof of Proposition 7.2. |
| Not addressed | Interlacing at every genus, the joint degree–genus regime, the b-deformed family, the multi-profile third term and the next-order occupancy correction are listed as open in Section 8. |

## Scope of assurance

Internal AI-assisted review only: written proofs checked by AI reviewers, exact computer checks replayed, sources read.
No human specialist review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
Literature searches were of limited scope (see [RESEARCH_GATES.json](RESEARCH_GATES.json)).
