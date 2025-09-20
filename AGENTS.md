# Thesis in Multivariate Conditional Density Estimation: Comparison with Transformation Forest Models, Copulas and Normalizing Flows

This repository collects the R notebooks and scripts used in the associated thesis.  All experiments are run directly in R to study multivariate conditional densities with mathematical rigour.

# Hardware environment

- macOS (MacBook Pro, Apple M1)
- 10 CPU cores available
- 16 GB RAM
- No discrete GPU acceleration

# Scientific approach
This repository is a research notebook rather than a software project. All scripts are executed directly in R to study problems in mathematical statistics and probability theory. The focus is on precise and reproducible experiments with **Mathematical Rigour**, not on standard software engineering workflows or packaging. The choice of R reflects the statistical setting; Python is intentionally avoided. Always analyze outputs for numerical validation. 

# Thesis Editing Rules for `report.rmd`


<thesis_agent>


  <scope>
    - Primary file: ./MSc_Report.Rnw  (rename from .rmd if needed)
    - Child chapters: chapter01.Rnw … chapter05.Rnw, then appendix files (chapterA*.Rnw)
    - Goal: clear scientific English, five chapters + Appendix, LaTeX/Sweave‑correct cross‑refs, build‑ready with knitr
  </scope>



  <editing_workflow>
    1. Read the whole file; map headings and cross-references.
    2. Plan edits; create a short checklist in comments near the changed section.
    3. Edit text first, then captions, then cross-refs. Keep math and code intact.
    4. After edits, scan for broken labels, undefined acronyms, and citation warnings.
  </editing_workflow>

<toc_five_chapters_plus_appendix>
    - Target ToC order:
      1) Introduction to the Research Problem
      2) Literature/Methodological Review or Background
      3) Data Analysis / Modeling / Validation / Simulation
      4) Interpretation of Results / Conclusion
      5) References
      Appendix (Code and materials)
  </toc_five_chapters_plus_appendix>

<cross_refs_and_labels>
    - Chapters: \chapter{...}\label{ch:intro|ch:background|ch:analysis|ch:interpretation|ch:refs}
    - Sections, figures, tables, algos, equations:
      * \label{sec:*}, \label{fig:*}, \label{tab:*}, \label{alg:*}, \label{eq:*}
      * Refer with \ref{fig:*}, \ref{tab:*}, \ref{alg:*}, and \eqref{eq:*}
    - Do not use R Markdown cross‑refs like (#fig:label)
    - Keep math in LaTeX; number equations only when referenced
  </cross_refs_and_labels>

  <citations>
    - Use natbib commands consistently: \citep for parenthetical, \citet for textual
    - No bare URLs in body text; put links in .bib
    - Ensure all cited keys exist in biblio.bib; no unused entries
  </citations>


  <global_style>
    <register_and_tone>
      - Formal, neutral, precise. No colloquialisms, metaphors, or humor.
      - American English for spelling and punctuation.
      - Use the Oxford comma; avoid exclamation marks; avoid italics for emphasis.
    </register_and_tone>

    <sentence_architecture>
      - One claim per sentence; target 12–20 words; split if >28 words.
      - Topic-before-new: start with known context, end with new information.
      - Prefer SVO order (Subject–Verb–Object). Avoid stacked subordinate clauses.
      - No rhetorical questions. Rewrite as declarative statements.
      - Place the main clause before qualifying clauses (“because…”, “although…”).
    </sentence_architecture>

    <voice_and_clarity>
      - Prefer active voice. Use passive only when the actor is unknown or irrelevant.
      - Replace nominalizations with verbs: “provide an analysis of” → “analyze”.
      - Replace vague intensifiers: “significantly”, “substantially”, “remarkably” → quantify or remove.
      - Replace deictics “this”, “that”, “it” with explicit nouns when ambiguous.
    </voice_and_clarity>

    <hedging_and_strength_of_claims>
      - Use calibrated verbs: “suggest”, “indicate”, “provide evidence”, “is consistent with”.
      - Reserve “prove”, “establish”, “demonstrate” for logically airtight arguments.
      - Use “significant” only for statistical significance and state the test/context if mentioned.
      - Prefer modal verbs to mark uncertainty: “may”, “can”, “is likely to”.
    </hedging_and_strength_of_claims>

    <cohesion_and_flow>
      - Each paragraph: Topic sentence → Evidence/Reasoning → Local takeaway/Bridge.
      - Transitional tokens allowed and limited: “However”, “In contrast”, “Moreover”, “Specifically”, “Consequently”.
      - Keep parallelism in lists and multi-part contributions.
    </cohesion_and_flow>

    <terminology_and_acronyms>
      - Define at first use: “Transformation Random Forests (TRF)”. Thereafter, use the acronym.
      - Maintain a single canonical term per concept. Add a glossary/notation table if needed.
      - Hyphenate compound modifiers before nouns: “lower‑triangular map”, “per‑sample time”.
    </terminology_and_acronyms>

    <numbers_units_and_symbols>
      - Use numerals for all measurements, metrics, counts, and statistics.
      - Use SI units and spaces between number and unit: “10 ms”, “3 GB”.
      - Keep consistent decimal precision across comparable values.
      - Spell out numbers only when they start a sentence; otherwise, rewrite the sentence.
    </numbers_units_and_symbols>

    <punctuation_and_formatting>
      - Commas: use for nonrestrictive clauses; avoid serial commas inside math.
      - Colons introduce lists or clarifications; semicolons join closely related clauses.
      - Avoid em dashes; prefer commas or parentheses sparingly.
      - Quotation marks only for exact quotes or defined terms at first introduction.
    </punctuation_and_formatting>

    <figure_table_reference_style>
      - Use “Figure 1”, “Table 2”, “Algorithm 3” consistently; no abbreviations.
      - Always cite a figure/table in running text before or near its appearance.
      - Captions summarize dataset, metric, and takeaway in 2–3 sentences.
    </figure_table_reference_style>

    <word_choice_macros>
      - Replace “very”, “really”, “quite”, “highly” → remove or quantify.
      - Replace “in order to” → “to”; “due to the fact that” → “because”.
      - Replace “different” when unspecified → name the concrete contrast.
      - Replace “respectively” unless mapping is unambiguous → rewrite explicitly.
      - Replace “since” for causality → “because”; keep “since” for time.
    </word_choice_macros>

    <consistency_checks_to_run>
      - Consistent US spelling (e.g., “modeling”, “behavior”, “labeled”).
      - Uniform hyphenation of repeated terms.
      - No mixed tenses within a paragraph unless justified by time scope.
      - No undefined symbols or acronyms.
      - No bare URLs in body text; use citation keys.
    </consistency_checks_to_run>

    <find_replace_snippets>
      - “There is|There are” + weak verb → rewrite with concrete subject.
      - “It is|It was [adjective] that …” → rewrite active with specific subject.
      - “This shows that …” after a number → include magnitude and unit.
      - Replace stacked parentheses with commas; keep at most one level of parentheses.
    </find_replace_snippets>
  </global_style>

  <tense_and_person>
    <person>
      - Use “we” or “this thesis”; do not use “I” or address the reader as “you”.
      - Use “this thesis” for document-structure statements; use “we” for author actions.
      - Do not use em dash or " — "
    </person>
 <didactic_patterns>
    - Introduction includes audience framing and a brief Non‑Goals paragraph
    - After key displayed equations, add a 2–3 sentence plain‑language paraphrase
    - For core operations, add a short numbered 1–3 step list
    - Early notation table in Chapter 2; define roles before depth
    - Use named bullet lists with bold lead‑ins for properties or reasons
    - Use “First/Then/Finally” for procedures
    - Name caveats and failure modes briefly; hedge, then give a takeaway
  </didactic_patterns>
    <tense_map_by_content>
      - Document facts and structure (what the text does): Present.
        * “This thesis investigates …” / “Section 3 presents …”
      - Definitions, theorems, equations, and general truths: Present.
        * “We define …”, “Equation (2) gives …”
      - Methods you designed or implemented: Past for what you did; Present for enduring properties.
        * “We implemented the training routine …” and “The algorithm runs in linear time.”
      - Experimental setup and execution: Past.
        * “We standardized features using training statistics only.”
      - Results observed in experiments: Past for measurements; Present for interpretation and figure/table references.
        * “Our model achieved a median test NLL of …” / “Table 2 shows the breakdown …”
      - Discussion of implications and limitations: Present with modals.
        * “The approach may degrade when …”, “This limitation likely stems from …”
      - Prior work (survey): Present for claims and established methods; Present perfect to summarize a research line.
        * “Smith (2020) proposes …” / “Prior studies have explored …”
      - Conclusion and contributions: Present; Future only for explicit future work.
        * “We provide …” / “Future work will examine …”
    </tense_map_by_content>

    <micro_rules>
      - Keep tense consistent within sentences and paragraphs unless time changes require a switch.
      - Refer to figures/tables in Present: “Figure 3 reports …”
      - Avoid “we show that …” in past; prefer Present: “we show that …”
      - Use Present perfect to connect past literature to current state: “has been shown”.
      - Do not mix Future with speculative language redundantly (“will likely” → “is likely to”).
    </micro_rules>

    <rewrite_examples>
      - Weak: “In this section, we will describe the method.” → Strong: “This section describes the method.”
      - Mixed tenses: “We evaluate and found improvements.” → “We evaluated and found improvements.”
      - Ambiguous: “This is significant.” → “The improvement is statistically significant (…);”
      - Passive to active: “Experiments were conducted on three datasets.” → “We conducted experiments on three datasets.”
    </rewrite_examples>
  </tense_and_person>


  <math_and_notation>
    - Define symbols before use; keep symbols and units consistent..
    - Do not start sentences with symbols; place punctuation outside math.
  </math_and_notation>

  <figures_and_tables>
    - Consistent model/order across all visuals; labeled axes and units.
    - Tables align decimals; include units and notes; avoid vertical rules.
    - Every figure/table is cited in text and has a self-contained caption.
  </figures_and_tables>

  <citations_and_references>
    - Cite primary sources for definitions; cite datasets and baselines at first mention.
    - Deduplicate references; ensure consistent capitalization in titles per chosen style.
  </citations_and_references>

  <language_quality_checks>
    - Articles present; concrete subject; precise verb.
    - Fix frequent German→English issues (articles, comma rules, decimal points).
    - Avoid mixed UK/US spelling; prefer US consistently.
  </language_quality_checks>

  <revision_policy>
    - Prefer surgical edits; do not change meaning or introduce new claims.
    - Do not add external citations or new results without explicit instruction.
    - Leave TODO comments only where information is missing and clearly mark them.
  </revision_policy>

  <commit_message_template>
    - refactor(text): tighten language in <section> (active voice, one claim per sentence)
    - fix(refs): repair labels/cross-refs and captions
    - chore(style): unify terminology, hyphenation, and acronyms
    - docs(results): clarify claims and add evidence links to existing tables/figures
  </commit_message_template>

  <quality_gate_checklist>
    - [ ] Abstract has Problem, Approach, Data, Results, Implication.
    - [ ] Each section starts with a topic sentence and ends with a bridge.
    - [ ] Active voice predominates; hedging used where appropriate.
    - [ ] All figures/tables referenced; captions self-contained.
    - [ ] Notation consistent; acronyms defined once.
    - [ ] Citations compile; references deduplicated.
    - [ ] Cross-refs resolve; knitting warnings addressed.
    - [ ] No R Markdown cross‑ref syntax remains
    - [ ] Figures/tables have captions and \label; all are referenced in text
    - [ ] All \ref and \eqref resolve; no LaTeX warnings for undefined labels
    - [ ] All \citep/\citet keys resolve; BibTeX runs clean
  </quality_gate_checklist>

</thesis_agent>

