"""
pyFroger — Atelier de groupes (interface Streamlit)
===================================================

Qualification des groupes de fautes alimentant la méthode par les groupes
(Froger). Chaque groupe variant est dépliable : on y lit *toutes* ses leçons
dans leur contexte, puis on déclare où porte la faute — sur la **variante**, ou
sur la **référence** (côté complémentaire, témoin de référence inclus). C'est
ce second choix qui permet de déplacer la racine hors de la référence : le
témoin *sans faute* devient l'archétype. La racine est donc déduite des fautes
déclarées, plus assignée d'office au texte de base.

Lancement :
    pip install streamlit matplotlib networkx pandas rich
    streamlit run app.py

L'appli importe ton vrai code pyfroger. La seule extension est la méthode
publique `UnorientedStemma.build_from_faults`, qui reconstruit le treillis à
partir des côtés fautifs choisis (la référence peut y figurer).
"""
from __future__ import annotations

import io
import csv
import tempfile
from collections import defaultdict

import pandas as pd
import streamlit as st

from pyfroger.list_apparatus import ListApparatus
from pyfroger.froger_viz import (PALETTE, sigla, build_stemma_from_faults,
                                 draw_stemma, draw_lattice)


st.set_page_config(page_title="pyFroger — Atelier de groupes",
                   page_icon="❧", layout="wide")

# ── Rôles possibles d'un groupe variant ──
ROLE_VAR = "la variante"
ROLE_REF = "la référence (complément)"
ROLE_POLY = "polygénèse"
ROLE_NEUT = "neutre"
ROLES_MULTI = [ROLE_VAR, ROLE_REF, ROLE_POLY, ROLE_NEUT]
ROLES_SINGULAR = [ROLE_VAR, ROLE_REF, ROLE_NEUT]
ROLE_SHORT = {ROLE_VAR: "faute · variante", ROLE_REF: "faute · référence",
              ROLE_POLY: "polygénèse", ROLE_NEUT: "neutre"}

st.markdown(
    f"""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=EB+Garamond:ital,wght@0,400;0,500;1,400&family=Inter:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');
    .stApp {{ background: {PALETTE['paper']}; }}
    html, body, [class*="css"] {{ font-family: 'Inter', sans-serif; color: {PALETTE['ink']}; }}
    .ed-eyebrow {{ font: 500 12px 'Inter'; letter-spacing: .18em; text-transform: uppercase;
                   color: {PALETTE['muted']}; margin-bottom: .2rem; }}
    .ed-title  {{ font: 500 2.4rem 'EB Garamond', serif; line-height: 1.05; margin: 0; }}
    .ed-sub    {{ font: italic 400 1.05rem 'EB Garamond', serif; color: {PALETTE['muted']};
                  margin-top: .35rem; }}
    .ed-rule   {{ border: none; border-top: 1px solid {PALETTE['rule']}; margin: 1rem 0 1.4rem; }}
    .sec-head  {{ font: 500 1.35rem 'EB Garamond', serif; margin: .2rem 0 .6rem; }}
    .verdict {{ border-radius: 10px; padding: .85rem 1.05rem; font-size: .95rem;
                border: 1px solid; margin-bottom: .4rem; }}
    .v-clean {{ background: #eef6f0; border-color: {PALETTE['laminar']}; color: #234c39; }}
    .v-warn  {{ background: #fbeeea; border-color: {PALETTE['rubric']}; color: #7a2616; }}
    .v-clean b, .v-warn b {{ font-weight: 600; }}
    .sigil {{ font: 500 13px 'IBM Plex Mono', monospace; }}
    .lectio {{ font-family: 'EB Garamond', serif; }}
    .ref-read {{ color: {PALETTE['muted']}; font-style: italic; }}
    code, .stCode, pre {{ font-family: 'IBM Plex Mono', monospace !important; }}
    section[data-testid="stSidebar"] {{ background: #f7f1e6; border-right: 1px solid {PALETTE['rule']}; }}
    .stButton button {{ border: 1px solid {PALETTE['indigo']}; color: {PALETTE['indigo']};
                        background: transparent; border-radius: 8px; font-weight: 500; }}
    .stButton button:hover {{ background: {PALETTE['indigo']}; color: white; }}
    </style>
    """,
    unsafe_allow_html=True,
)


# ───────────────────────────── helpers (cache) ──────────────────────────────
@st.cache_data(show_spinner=False)
def read_header(csv_bytes: bytes, delimiter: str) -> list[str]:
    """Retourne les noms de colonnes (en-tête) du CSV de collation."""
    text = csv_bytes.decode("utf-8", errors="replace")
    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    try:
        return [h.strip() for h in next(reader)]
    except StopIteration:
        return []


@st.cache_data(show_spinner=False)
def load_apparatus(csv_bytes: bytes, ref: str, delimiter: str):
    with tempfile.NamedTemporaryFile("wb", suffix=".csv", delete=False) as fh:
        fh.write(csv_bytes)
        path = fh.name
    return ListApparatus.from_csv(path, ref_witness=ref, delimiter=delimiter)


def build_locus_context(app: ListApparatus):
    """Reconstruit, pour chaque lieu variant, la distribution des leçons.

    Renvoie deux dictionnaires :
        dist : locus -> {leçon -> liste triée de témoins}
        refr : locus -> leçon du témoin de référence
    """
    dist: dict[str, dict[str, list[str]]] = defaultdict(dict)
    refr: dict[str, str] = {}
    for g in app.groups:
        dist[g.locus][g.reading] = sorted(g.witnesses)
        refr[g.locus] = g.reference
    return dist, refr


def context_line(locus: str, dist, refr, ref_name: str) -> str:
    """Une ligne lisible : toutes les leçons d'un lieu, réf. indiquée."""
    parts = []
    for reading, wits in sorted(dist.get(locus, {}).items(),
                                key=lambda kv: (-len(kv[1]), kv[0])):
        parts.append(f"« {reading} » ({', '.join(wits)})")
    ref_read = refr.get(locus, "")
    body = "   ·   ".join(parts)
    return f"{body}    [réf. {ref_name} : « {ref_read} »]"


# ───────────────────────────────── header ───────────────────────────────────
st.markdown(
    """
    <h1 class="ed-title">Méthode par les groupes</h1>
    <div class="ed-sub">Qualification des fautes communes &amp; orientation
    déduite : la racine n'est plus assignée d'office à la référence.</div>
    <hr class="ed-rule"/>
    """,
    unsafe_allow_html=True,
)


# ──────────────────────────────── Barre latérale ────────────────────────────
with st.sidebar:
    st.markdown('<div class="sec-head">Collation</div>',
                unsafe_allow_html=True)
    source = st.radio("Source", ["Exemple artificiel", "Importer un CSV"],
                      label_visibility="collapsed")

    if source.startswith("Importer"):
        up = st.file_uploader("Fichier de collation (.csv)", type=["csv"])
        delim = st.text_input("Délimiteur", value=",", max_chars=1) or ","
        if up is None:
            st.info("Importe un CSV pour commencer.")
            st.stop()
        csv_bytes = up.read()
    else:
        delim = ","
        with open("demo_collation.csv", "rb") as fh:
            csv_bytes = fh.read()

    columns = read_header(csv_bytes, delim)
    witness_cols = columns[1:] if len(columns) > 1 else []
    if not witness_cols:
        st.error("En-tête illisible : aucune colonne de témoin trouvée.")
        st.stop()

    default_ref = "T1" if "T1" in witness_cols else witness_cols[0]

    st.markdown("<hr class='ed-rule'/>", unsafe_allow_html=True)
    st.markdown('<div class="sec-head">Texte de base / référence</div>',
                unsafe_allow_html=True)
    st.caption("Témoin contre lequel les leçons sont collationnées. Il définit "
               "les *variantes*, mais n'est plus forcément la racine de l'arbre "
               "— celle-ci se déduit des fautes déclarées.")
    ref = st.selectbox("Témoin de collation", witness_cols,
                       index=witness_cols.index(default_ref),
                       label_visibility="collapsed")

    try:
        app = load_apparatus(csv_bytes, ref, delim)
    except Exception as exc:                       # noqa: BLE001
        st.error(f"Lecture impossible : {exc}")
        st.stop()

    full_witnesses = set(app.witnesses)
    universe = frozenset(full_witnesses | {ref})
    freq = app.group_by_frequency()
    if not freq:
        st.error("Aucun groupe de fautes trouvé (toutes les leçons égalent la "
                 "référence ?).")
        st.stop()

    dist, refr = build_locus_context(app)
    all_keys = list(freq.keys())

    def complement(k):
        return frozenset(universe - k)

    n_loci_total = sum(v["frequency"] for v in freq.values())
    st.markdown(
        f'<div style="font-size:.85rem;color:{PALETTE["muted"]};margin-top:.4rem">'
        f"Témoins&nbsp;: <span class='sigil'>{', '.join(sorted(full_witnesses))}</span>"
        f" &nbsp;·&nbsp; base&nbsp;: <span class='sigil'>{ref}</span><br>"
        f"{len(freq)} groupes variants · {n_loci_total} lieux de variants</div>",
        unsafe_allow_html=True,
    )

    # ── État de qualification : un rôle par groupe (source de vérité) ──
    fingerprint = (ref, tuple(sorted(sigla(k) for k in all_keys)))
    if st.session_state.get("loaded_ref") != fingerprint:
        st.session_state.role = {k: ROLE_VAR for k in all_keys}
        st.session_state.loaded_ref = fingerprint
        st.session_state.nonce = 0
    role = st.session_state.role

    def side_of(k):
        """Côté fautif (frozenset) choisi pour le groupe k, ou None si écarté."""
        r = role.get(k, ROLE_VAR)
        if r == ROLE_VAR:
            return k
        if r == ROLE_REF:
            return complement(k)
        return None

    def current_sides():
        return [s for s in (side_of(k) for k in all_keys) if s is not None]

    def side_to_key():
        m = {}
        for k in all_keys:
            s = side_of(k)
            if s is not None:
                m[s] = k
        return m

    def overlapping_group_keys():
        try:
            tmp = build_stemma_from_faults(app, current_sides(),
                                           full_witnesses, ref)
        except ValueError:
            return set()
        s2k = side_to_key()
        keys = set()
        for g1, g2 in (tmp.contamination_overlap or []):
            for node in (g1, g2):
                kk = s2k.get(node)
                if kk is not None:
                    keys.add(kk)
        return keys

    st.markdown("<hr class='ed-rule'/>", unsafe_allow_html=True)
    st.markdown('<div class="sec-head">Présélection</div>',
                unsafe_allow_html=True)
    st.caption("Réglages en lot ; affine ensuite groupe par groupe ci-contre.")

    def set_roles(fn):
        for k in all_keys:
            role[k] = fn(k)
        st.session_state.nonce += 1

    c1, c2 = st.columns(2)
    with c1:
        if st.button("Tout = variante", use_container_width=True):
            set_roles(lambda k: ROLE_VAR)
            st.rerun()
        if st.button("Garder ≥2", use_container_width=True,
                     help="Mettre les leçons singulières en « neutre »."):
            set_roles(lambda k: ROLE_VAR if len(k) >= 2 else ROLE_NEUT)
            st.rerun()
    with c2:
        if st.button("Sans conflit", use_container_width=True):
            confl = overlapping_group_keys()
            set_roles(lambda k: ROLE_NEUT if k in confl else role.get(k, ROLE_VAR))
            st.rerun()
        if st.button("Réinitialiser", use_container_width=True):
            set_roles(lambda k: ROLE_VAR)
            st.rerun()

    if st.button("Conflits → polygénèse", use_container_width=True):
        confl = overlapping_group_keys()
        set_roles(lambda k: ROLE_POLY if (k in confl and len(k) >= 2)
                  else role.get(k, ROLE_VAR))
        st.rerun()
    st.caption("Marque les accords qui se chevauchent comme coïncidences "
               "(hasard) : ils sortent de la construction.")

    max_freq = max((v["frequency"] for v in freq.values()), default=1)
    min_loci = st.slider(
        "Seuil « fort » (loci)", 1, max(2, max_freq), min(3, max(2, max_freq)),
        help="En deçà, un accord est jugé « polygénèse plausible / faible ».")

    with st.expander("Avancé — forcer la racine"):
        forced = st.selectbox(
            "Racine", ["(déduite des fautes)"] + sorted(universe),
            help="Normalement la racine est déduite (témoin sans faute). "
                 "Forcer une racine peut produire un arbre incohérent si ce "
                 "témoin porte des fautes.")
        forced_root = None if forced.startswith("(") else forced


# ───────────────────────────────── Corps ────────────────────────────────────
left, right = st.columns([0.46, 0.54], gap="large")

with left:
    st.markdown('<div class="sec-head">Groupes variants — qualification</div>',
                unsafe_allow_html=True)
    st.caption("Déplie un groupe pour lire toutes ses leçons en contexte, puis "
               "déclare **où porte la faute** : sur *la variante* (les témoins "
               "du groupe ont innové) ou sur *la référence* (le côté "
               "complémentaire, référence comprise, a innové — ce qui éloigne "
               "la racine de la référence). **Polygénèse** / **neutre** ⇒ "
               "écarté. Le nombre d'accords (en-tête) mesure la force de "
               "l'indice : un seul accord se dissout volontiers en hasard.")

    hide_singulars = st.checkbox("Masquer les leçons singulières (|G|=1)",
                                 value=False)

    nonce = st.session_state.nonce
    order = sorted(all_keys,
                   key=lambda s: (-len(s), -freq[s]["frequency"], sorted(s)))

    for k in order:
        if hide_singulars and len(k) == 1:
            continue
        s = sigla(k)
        info = freq[k]
        n = info["frequency"]
        singular = len(k) == 1
        comp = complement(k)

        r = role.get(k, ROLE_VAR)
        excluded = r in (ROLE_POLY, ROLE_NEUT)
        mark = "∅ " if excluded else ("⇲ " if r == ROLE_REF else "")
        head = (f"{mark}{s}   ·   |G|={len(k)}   ·   {n} accords c. réf.   "
                f"·   [{ROLE_SHORT.get(r, r)}]")

        with st.expander(head, expanded=False):
            strength = ("base solide" if n >= min_loci
                        else "base faible — polygénèse plausible")
            st.caption(f"{n} lieu(x) d'accord — {strength}.")
            st.markdown(
                f"<div style='font-size:.85rem'>"
                f"variante : <span class='sigil'>{s}</span>"
                f" &nbsp;·&nbsp; référence (complément) : "
                f"<span class='sigil'>{sigla(comp)}</span></div>",
                unsafe_allow_html=True)

            opts = ROLES_SINGULAR if singular else ROLES_MULTI
            cur = r if r in opts else ROLE_VAR
            chosen = st.selectbox(
                "La faute porte sur", opts, index=opts.index(cur),
                key=f"role_{s}_{nonce}")
            role[k] = chosen
            if singular:
                st.caption("Leçon singulière : « variante » = ce témoin a innové "
                           "(simple feuille) ; « référence » = lui seul conserve "
                           "l'original — il devient très ancestral.")

            detail = []
            for locus, reading, refread in zip(
                    info["loci"], info["readings"], info["reference"]):
                detail.append({
                    "locus": locus,
                    "leçon du groupe": reading,
                    "référence": refread,
                    "contexte (toutes leçons)": context_line(
                        locus, dist, refr, ref),
                })
            st.dataframe(pd.DataFrame(detail), hide_index=True,
                         use_container_width=True,
                         column_config={
                             "contexte (toutes leçons)":
                                 st.column_config.TextColumn(width="large"),
                         })

# Les rôles déclarés = paramètres de construction.
selected_keys = {k for k in all_keys if side_of(k) is not None}
faulty_sides = current_sides()
s2k = side_to_key()

with right:
    if not selected_keys:
        st.info("Déclare au moins une faute (variante ou référence) pour "
                "construire le stemma.")
        st.stop()
    try:
        stemma = build_stemma_from_faults(app, faulty_sides, full_witnesses,
                                          ref, root=forced_root)
    except ValueError as exc:
        st.markdown(
            f'<div class="verdict v-warn">▲ <b>Fautes contradictoires.</b> '
            f"{exc}</div>", unsafe_allow_html=True)
        st.stop()
    except Exception as exc:                       # noqa: BLE001
        st.error(f"Construction impossible avec cette sélection : {exc}")
        st.stop()

    overlaps = stemma.contamination_overlap or []
    contaminated = sorted(stemma.contaminated_witnesses or [])
    root_found = getattr(stemma, "root", ref)
    root_kind = "forcée" if forced_root else "déduite"

    if not overlaps:
        st.markdown(
            f'<div class="verdict v-clean">✦ <b>Stemma laminaire.</b> '
            f"Racine {root_kind} : <span class='sigil'>{root_found}</span>. "
            "Tous les groupes s'emboîtent — transmission verticale pure.</div>",
            unsafe_allow_html=True)
    else:
        st.markdown(
            f'<div class="verdict v-warn">▲ <b>{len(overlaps)} conflit(s) '
            f"d'emboîtement.</b> Racine {root_kind} : "
            f"<span class='sigil'>{root_found}</span>. Témoins impliqués : "
            f"{', '.join(contaminated)}.</div>", unsafe_allow_html=True)

    m1, m2, m3 = st.columns(3)
    m1.metric("Fautes retenues", len(selected_keys))
    m2.metric("Conflits", len(overlaps))
    m3.metric("Témoins contaminés", len(contaminated))

    tab_tree, tab_lat, tab_conf, tab_nwk = st.tabs(
        ["Stemma", "Emboîtement", "Conflits", "Newick"])

    with tab_tree:
        st.caption(f"Arbre enraciné sur « {root_found} » ({root_kind} des "
                   "fautes déclarées).")
        st.pyplot(draw_stemma(stemma.tree), use_container_width=True)
        with st.expander("Arbre ASCII (sortie native de pyfroger)"):
            st.code(stemma.display_tree(), language="text")

    with tab_lat:
        st.caption("Empilement par inclusivité. Indigo = couverture. "
                   "Rouge pointillé = chevauchement non emboîté (témoin partagé "
                   "annoté).")
        st.pyplot(draw_lattice(stemma), use_container_width=True)

    # ── Conflits : variantes à arbitrer (features 5 + 6) ──
    with tab_conf:
        if not overlaps:
            st.success("Aucun conflit : rien à arbitrer.")
        else:
            st.caption("Chaque côté fautif en conflit est détaillé lieu par "
                       "lieu. Tague-le **polygénèse** (hasard) pour l'écarter, "
                       "ou garde-le comme **contamination réelle**. "
                       f"Force « forte » si loci ≥ {min_loci}.")
            seen = set()
            for g1, g2 in overlaps:
                shared = ", ".join(sorted(g1 & g2))
                st.markdown(
                    f"<div style='margin:.3rem 0 .1rem'>"
                    f"<span class='sigil'>{sigla(g1)}</span> ⋈ "
                    f"<span class='sigil'>{sigla(g2)}</span> "
                    f"<span style='color:{PALETTE['rubric']}'>"
                    f"(témoin(s) partagé(s) : {shared})</span></div>",
                    unsafe_allow_html=True)
                for node in (g1, g2):
                    k = s2k.get(node)
                    if k is None or k in seen:
                        continue
                    seen.add(k)
                    info = freq.get(k, {})
                    n = info.get("frequency", 0)
                    force = ("forte" if n >= min_loci
                             else "faible · polygénèse possible")
                    side_lbl = ("référence" if role.get(k) == ROLE_REF
                                else "variante")
                    cols = st.columns([0.62, 0.18, 0.20])
                    with cols[0]:
                        st.markdown(
                            f"**{sigla(k)}** · côté *{side_lbl}* — {n} loci · "
                            f"*{force}*")
                        ev = []
                        for locus, reading, refread in zip(
                                info.get("loci", []),
                                info.get("readings", []),
                                info.get("reference", [])):
                            ev.append(f"`{locus}` « {reading} » ← « {refread} »")
                        st.markdown(" · ".join(ev) if ev else "—")
                    with cols[1]:
                        if st.button("→ polygénèse", key=f"cpoly_{sigla(k)}",
                                     use_container_width=True):
                            role[k] = ROLE_POLY
                            st.session_state.nonce += 1
                            st.rerun()
                    with cols[2]:
                        if st.button("contamination", key=f"ckeep_{sigla(k)}",
                                     use_container_width=True,
                                     help="Garder : chevauchement réel "
                                          "(transmission horizontale)."):
                            st.session_state.nonce += 1
                            st.rerun()
                st.markdown("<hr class='ed-rule'/>", unsafe_allow_html=True)
            st.caption("Astuce : les accords reposant sur peu de loci sont les "
                       "premiers candidats à la polygénèse.")

    with tab_nwk:
        st.code(stemma.to_newick(), language="text")
        st.caption(f"Sérialisation enracinée sur « {root_found} » "
                   f"({root_kind} des fautes déclarées).")


st.markdown("<hr class='ed-rule'/>", unsafe_allow_html=True)
with st.expander("Méthode — où porte la faute, et pourquoi la racine bouge"):
    st.markdown(
        "Une variation oppose deux camps : les porteurs de la *variante* et "
        "ceux qui lisent avec la *référence*. Dans un arbre **non orienté**, "
        "rien ne dit a priori lequel est l'innovation : le déclarer, c'est "
        "orienter. Tant que la faute est toujours mise « sur la variante », la "
        "référence ne porte jamais de faute et reste mécaniquement la racine — "
        "la *référence imaginaire*. En autorisant la faute « sur la référence "
        "» (le **complément**, référence comprise), ce témoin peut être éliminé "
        "à son tour : la racine devient le témoin réellement sans faute. "
        "Techniquement, `build_from_faults()` reconstruit le treillis de "
        "couverture à partir des seuls côtés fautifs retenus, sans présupposer "
        "que la référence est ancestrale. Seules les **fautes communes** "
        "fondent un groupe ; le nombre de **loci** partagés jauge le risque de "
        "polygénèse (un accord isolé se dissout volontiers en hasard, une "
        "dizaine non).")