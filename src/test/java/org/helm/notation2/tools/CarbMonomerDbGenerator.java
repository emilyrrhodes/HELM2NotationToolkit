package org.helm.notation2.tools;

import java.io.FileOutputStream;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import org.helm.chemtoolkit.AbstractChemistryManipulator.StType;
import org.helm.notation2.Attachment;
import org.helm.notation2.Chemistry;
import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.calculation.MoleculePropertyCalculator;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.tools.MonomerParser;
import org.helm.notation2.wsadapter.MonomerStoreConfiguration;
import org.jdom2.Attribute;
import org.jdom2.Element;
import org.jdom2.output.Format;
import org.jdom2.output.XMLOutputter;

/**
 * One-off, reproducible generator for the CARB section of the shipped monomer
 * database ({@code MonomerDBGZEncoded.xml}).
 *
 * <p>Each monosaccharide is defined here as an R-group-annotated canonical
 * SMILES following the glycosidic-bond cap convention already verified by
 * {@link CarbBuilderMoleculeTest} (R1/anomeric: the whole hydroxyl is the cap,
 * "O[*:n]", exposing the bare ring carbon; every other free hydroxyl: only the
 * hydroxyl hydrogen is the cap, "[*:n][H]", so the bridging oxygen is retained;
 * net loss per glycosidic bond is one H2O). Each stereochemistry was verified,
 * before being committed here, by stripping the atom-map annotations and
 * confirming the resulting neutral structure canonicalises identically to the
 * authoritative PubChem isomeric SMILES for that exact anomer/configuration.
 *
 * <p>Running {@link #main} regenerates {@code target/carb-polymer-block.xml},
 * whose {@code <Polymer polymerType="CARB">} element is copied into the shipped
 * resource. It also self-checks every monomer: the molfile is generated from the
 * SMILES by the chemistry engine, the monomer passes {@link MonomerParser#validateMonomer},
 * and a single free monomer builds to its known molecular formula.
 *
 * <p>This class is intentionally not named {@code *Test}, so it is skipped by the
 * normal Surefire run. Regenerate with, from the module root:
 * <pre>
 *   mvn -o test-compile
 *   java -cp "target/test-classes:target/classes:$(cat /tmp/cp.txt)" \
 *        org.helm.notation2.tools.CarbMonomerDbGenerator
 * </pre>
 */
public final class CarbMonomerDbGenerator {

  /** OH cap on the anomeric carbon: the whole hydroxyl leaves on bond formation. */
  private static final String OH = Attachment.CAP_GROUP_OH;
  /** H cap on a ring/primary hydroxyl: only the H leaves, the O bridges. */
  private static final String H = Attachment.CAP_GROUP_H;

  private CarbMonomerDbGenerator() {
  }

  /** A single monosaccharide definition and its expected free-monomer formula. */
  private static final class Sugar {
    final String id;
    final String name;
    final String smiles;
    final String formula;
    /** attachment labels, in order; the first entry is the anomeric (OH) cap. */
    final String[] labels;

    Sugar(String id, String name, String smiles, String formula, String... labels) {
      this.id = id;
      this.name = name;
      this.smiles = smiles;
      this.formula = formula;
      this.labels = labels;
    }
  }

  private static List<Sugar> sugars() {
    List<Sugar> s = new ArrayList<Sugar>();
    // Hexoses - R1 anomeric, R2/R3/R4 ring hydroxyls, R6 primary alcohol.
    s.add(new Sugar("b-D-Glcp", "beta-D-Glucopyranose",
        "[H:6]OC[C@H]1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    s.add(new Sugar("a-D-Glcp", "alpha-D-Glucopyranose",
        "[H:6]OC[C@H]1O[C@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    s.add(new Sugar("b-D-Galp", "beta-D-Galactopyranose",
        "[H:6]OC[C@H]1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    s.add(new Sugar("a-D-Galp", "alpha-D-Galactopyranose",
        "[H:6]OC[C@H]1O[C@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    s.add(new Sugar("a-D-Manp", "alpha-D-Mannopyranose",
        "[H:6]OC[C@H]1O[C@H]([OH:1])[C@@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    s.add(new Sugar("b-D-Manp", "beta-D-Mannopyranose",
        "[H:6]OC[C@H]1O[C@@H]([OH:1])[C@@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H12O6", "R1", "R2", "R3", "R4", "R6"));
    // N-acetylhexosamines - C2 bears N-acetyl (no free hydroxyl, so no R2).
    s.add(new Sugar("b-D-GlcpNAc", "2-acetamido-2-deoxy-beta-D-glucopyranose",
        "[H:6]OC[C@H]1O[C@@H]([OH:1])[C@H](NC(C)=O)[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C8H15NO6", "R1", "R3", "R4", "R6"));
    s.add(new Sugar("a-D-GalpNAc", "2-acetamido-2-deoxy-alpha-D-galactopyranose",
        "[H:6]OC[C@H]1O[C@H]([OH:1])[C@H](NC(C)=O)[C@@H](O[H:3])[C@H]1O[H:4]",
        "C8H15NO6", "R1", "R3", "R4", "R6"));
    s.add(new Sugar("b-D-GalpNAc", "2-acetamido-2-deoxy-beta-D-galactopyranose",
        "[H:6]OC[C@H]1O[C@@H]([OH:1])[C@H](NC(C)=O)[C@@H](O[H:3])[C@H]1O[H:4]",
        "C8H15NO6", "R1", "R3", "R4", "R6"));
    // 6-deoxy-L-galactose (fucose) - C6 is a methyl, so no R6.
    s.add(new Sugar("a-L-Fucp", "alpha-L-Fucopyranose",
        "C[C@@H]1O[C@@H]([OH:1])[C@@H](O[H:2])[C@H](O[H:3])[C@@H]1O[H:4]",
        "C6H12O5", "R1", "R2", "R3", "R4"));
    // Pentose (xylose) - no C6 at all.
    s.add(new Sugar("b-D-Xylp", "beta-D-Xylopyranose",
        "C1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C5H10O5", "R1", "R2", "R3", "R4"));
    // Uronic acids - C6 is a carboxyl, so no R6.
    s.add(new Sugar("b-D-GlcpA", "beta-D-Glucopyranuronic acid",
        "OC(=O)[C@H]1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H10O7", "R1", "R2", "R3", "R4"));
    s.add(new Sugar("a-L-IdopA", "alpha-L-Idopyranuronic acid",
        "OC(=O)[C@@H]1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]",
        "C6H10O7", "R1", "R2", "R3", "R4"));
    // Sialic acid (Neu5Ac) - a ketose: the anomeric carbon is C2 (OH cap on R2),
    // C1 is a carboxyl, C3 is a deoxy methylene, C5 bears N-acetyl; free hydroxyls
    // at C4/C7/C8/C9.
    s.add(new Sugar("a-D-Neup5Ac", "alpha-N-acetylneuraminic acid",
        "CC(=O)N[C@@H]1[C@H](C[C@@](O[C@H]1[C@@H]([C@@H](CO[H:9])O[H:8])O[H:7])(C(=O)O)[OH:2])O[H:4]",
        "C11H19NO9", "R2", "R4", "R7", "R8", "R9"));
    return s;
  }

  public static void main(String[] args) throws Exception {
    setUseWebservice(false);

    // Build with NO explicit namespace, mirroring MonomerFactory#buildMonomerDbXMLFromCache
    // and MonomerParser#getMonomerElement: the shipped resource declares xmlns="lmr"
    // once on <MonomerDB>, and every element inherits it. Adding an explicit namespace
    // here would emit stray xmlns="" attributes and break loading.
    Element polymer = new Element("Polymer");
    polymer.setAttribute(new Attribute("polymerType", Monomer.CARBOHYDRATE_POLYMER_TYPE));

    for (Sugar sugar : sugars()) {
      Monomer m = buildMonomer(sugar);

      // 1) structural self-check: the monomer must satisfy the same validation
      //    the loader applies when reading the shipped resource.
      if (!MonomerParser.validateMonomer(m)) {
        throw new IllegalStateException("validateMonomer failed for " + sugar.id);
      }

      // 2) register and confirm a single free monomer builds to its known formula.
      registerCarb(m);
      String formula = MoleculePropertyCalculator.getMolecularFormular(
          HELM2NotationUtils.readNotation("CARB1{[" + sugar.id + "]}$$$$V2.0"));
      if (!sugar.formula.equals(formula)) {
        throw new IllegalStateException(
            sugar.id + " built to " + formula + ", expected " + sugar.formula);
      }
      System.out.println("OK " + sugar.id + " -> " + formula);

      polymer.getChildren().add(MonomerParser.getMonomerElement(m));
    }

    String out = "target/carb-polymer-block.xml";
    try (FileOutputStream fos = new FileOutputStream(out)) {
      // output just the <Polymer> subtree, to be spliced under <PolymerList> in
      // the shipped resource (which already provides the lmr default namespace).
      new XMLOutputter(Format.getPrettyFormat()).output(polymer, fos);
    }
    System.out.println("Wrote " + out);
  }

  private static Monomer buildMonomer(Sugar sugar) throws Exception {
    Monomer m = new Monomer();
    m.setAlternateId(sugar.id);
    m.setPolymerType(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    m.setMonomerType(Monomer.BACKBONE_MOMONER_TYPE);
    m.setName(sugar.name);
    m.setCanSMILES(sugar.smiles);
    // molfile is required to be present/valid but is not used by the CARB build
    // path (which builds from canSMILES); generate it from the SMILES. convert()'s
    // StType names the INPUT type, so StType.SMILES means "SMILES in, molfile out".
    m.setMolfile(Chemistry.getInstance().getManipulator().convert(sugar.smiles, StType.SMILES));

    for (int i = 0; i < sugar.labels.length; i++) {
      String label = sugar.labels[i];
      // the first label is the anomeric position (whole-hydroxyl OH cap); the
      // rest are ring/primary hydroxyls (H cap).
      boolean anomeric = (i == 0);
      String capName = anomeric ? OH : H;
      int n = Integer.parseInt(label.substring(1));
      String capSmiles = anomeric ? ("O[*:" + n + "]") : ("[*:" + n + "][H]");
      m.getAttachmentList().add(attachment(label, capName, capSmiles));
    }
    return m;
  }

  private static Attachment attachment(String label, String capGroupName, String capGroupSMILES) {
    Attachment att = new Attachment(label, capGroupName);
    att.setAlternateId(label + "-" + capGroupName);
    att.setCapGroupSMILES(capGroupSMILES);
    return att;
  }

  // addMonomer() silently skips an id already present, so evict any pre-existing
  // (possibly stale, cache-persisted) entry first before registering.
  private static void registerCarb(Monomer monomer) throws Exception {
    java.util.Map<String, java.util.Map<String, Monomer>> db =
        MonomerFactory.getInstance().getMonomerStore().getMonomerDB();
    java.util.Map<String, Monomer> carbMap = db.get(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    if (carbMap != null) {
      carbMap.remove(monomer.getAlternateId());
    }
    MonomerFactory.getInstance().getMonomerStore().addMonomer(monomer);
  }

  private static void setUseWebservice(boolean value) throws Exception {
    Field wsField = MonomerStoreConfiguration.class.getDeclaredField("isUseWebservice");
    wsField.setAccessible(true);
    wsField.set(MonomerStoreConfiguration.getInstance(), value);
  }
}
