package org.helm.notation2.tools;

import java.io.IOException;
import java.lang.reflect.Field;

import org.helm.chemtoolkit.CTKException;
import org.helm.notation2.Attachment;
import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.calculation.MoleculePropertyCalculator;
import org.helm.notation2.exception.BuilderMoleculeException;
import org.helm.notation2.exception.ChemistryException;
import org.helm.notation2.exception.ConnectionNotationException;
import org.helm.notation2.exception.GroupingNotationException;
import org.helm.notation2.exception.MonomerException;
import org.helm.notation2.exception.MonomerLoadingException;
import org.helm.notation2.exception.NotationException;
import org.helm.notation2.exception.ParserException;
import org.helm.notation2.exception.PolymerIDsException;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.wsadapter.MonomerStoreConfiguration;
import org.jdom2.JDOMException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * End-to-end structure-generation tests for CARB, using a real
 * D-glucopyranose ("Glcp") monomer definition to prove the full pipeline -
 * parsing, R-group validation, and {@code BuilderMolecule.buildMoleculefromCARB}
 * - produces chemically correct molecules, not just that the notation parses.
 *
 * The glycosidic-bond cap convention this monomer relies on (R1/anomeric:
 * whole hydroxyl is the cap, exposing the bare ring carbon; R2/R3/R4/R6:
 * only the hydroxyl's hydrogen is the cap, so the ring oxygen is retained and
 * becomes the bridging atom) was verified against the real "MCC" (OH cap)
 * and "hxy" (H cap) monomer entries already in MonomerDBGZEncoded.xml before
 * being applied here.
 */
public class CarbBuilderMoleculeTest {

  @BeforeClass
  public void registerGlucose() throws IOException, MonomerException, ChemistryException,
      NoSuchFieldException, IllegalAccessException {
    Field wsField = MonomerStoreConfiguration.class.getDeclaredField("isUseWebservice");
    wsField.setAccessible(true);
    wsField.set(MonomerStoreConfiguration.getInstance(), false);
    registerCarb(glucopyranose());
  }

  // addMonomer() silently skips an id already in the store, so evict any
  // pre-existing (possibly stale, cache-persisted) CARB entry first - otherwise a
  // leftover "b-D-Glcp" without these SMILES/cap definitions would make the build
  // tests fail depending on the shared ~/.helm monomer cache state.
  private void registerCarb(Monomer monomer) throws IOException, MonomerException, ChemistryException {
    java.util.Map<String, java.util.Map<String, Monomer>> db =
        MonomerFactory.getInstance().getMonomerStore().getMonomerDB();
    java.util.Map<String, Monomer> carbMap = db.get(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    if (carbMap != null) {
      carbMap.remove(monomer.getAlternateId());
    }
    MonomerFactory.getInstance().getMonomerStore().addMonomer(monomer);
  }

  @Test
  public void testSingleFreeMonomerHasGlucoseFormula() throws ParserException, JDOMException,
      BuilderMoleculeException, CTKException, NotationException, ChemistryException {
    // a single, fully unsubstituted CARB monomer should build to plain glucose
    String notation = "CARB1{[b-D-Glcp]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    String formula = MoleculePropertyCalculator.getMolecularFormular(helm2notation);
    Assert.assertEquals(formula, "C6H12O6");
  }

  @Test
  public void testUnknownLinkageIsRejectedAtBuildTime() throws ParserException, JDOMException,
      NotationException, ChemistryException, CTKException, PolymerIDsException, MonomerException,
      GroupingNotationException, ConnectionNotationException, MonomerLoadingException,
      org.helm.notation2.parser.exceptionparser.NotationException {
    // "R?" is valid notation for an unresolved/unknown linkage - it must pass
    // validation, but building a concrete molecule structure requires a known
    // attachment point, so build must fail with a clear, dedicated message.
    String notation = "CARB1{[b-D-Glcp].R?:[b-D-Glcp]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    try {
      MoleculePropertyCalculator.getMolecularFormular(helm2notation);
      Assert.fail("Expected BuilderMoleculeException for an unknown (R?) linkage");
    } catch (BuilderMoleculeException e) {
      Assert.assertTrue(e.getMessage().contains("unknown (R?) linkage"), e.getMessage());
    }
  }

  @Test
  public void testDisaccharideLosesOneWater() throws ParserException, JDOMException,
      BuilderMoleculeException, CTKException, NotationException, ChemistryException,
      PolymerIDsException, MonomerException, GroupingNotationException, ConnectionNotationException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // two glucose units joined by one glycosidic bond (1->4) should be a
    // standard disaccharide formula: 2x C6H12O6 minus one H2O = C12H22O11
    String notation = "CARB1{[b-D-Glcp].R4:[b-D-Glcp]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    String formula = MoleculePropertyCalculator.getMolecularFormular(helm2notation);
    Assert.assertEquals(formula, "C12H22O11");
  }

  @Test
  public void testBranchedGlycanLosesTwoWaters() throws ParserException, JDOMException,
      BuilderMoleculeException, CTKException, NotationException, ChemistryException,
      PolymerIDsException, MonomerException, GroupingNotationException, ConnectionNotationException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // A branched trisaccharide: a trunk of two glucoses joined (1->4), with a
    // third glucose branching onto the second monomer's R2. This exercises the
    // union-find branch-convergence path in BuilderMolecule.buildMoleculefromCARB
    // (the branch monomer merges back onto a trunk monomer that already carries an
    // incoming main-chain bond). Two glycosidic bonds: 3x C6H12O6 minus 2x H2O.
    String notation = "CARB1{[b-D-Glcp].R4:[b-D-Glcp]([b-D-Glcp].R2)}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    String formula = MoleculePropertyCalculator.getMolecularFormular(helm2notation);
    Assert.assertEquals(formula, "C18H32O16");
  }

  @Test
  public void testIntegerRepeatBuildsExpandedChain() throws ParserException, JDOMException,
      BuilderMoleculeException, CTKException, NotationException, ChemistryException,
      PolymerIDsException, MonomerException, GroupingNotationException, ConnectionNotationException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // An integer repeat "(R4:[b-D-Glcp])'2'" is expanded by the parser's numbering
    // into two concrete extra positions, each joined (1->4) to its predecessor. The
    // leading trunk monomer plus two repeat instances is a linear trisaccharide:
    // 3x C6H12O6 minus 2x H2O = C18H32O16.
    String notation = "CARB1{[b-D-Glcp].(R4:[b-D-Glcp])'2'}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    String formula = MoleculePropertyCalculator.getMolecularFormular(helm2notation);
    Assert.assertEquals(formula, "C18H32O16");
  }

  @Test
  public void testDirectMonomerRepeatBuildsExpandedChain() throws ParserException, JDOMException,
      BuilderMoleculeException, CTKException, NotationException, ChemistryException,
      PolymerIDsException, MonomerException, GroupingNotationException, ConnectionNotationException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    String notation = "CARB1{[b-D-Glcp]'3'}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    String formula = MoleculePropertyCalculator.getMolecularFormular(helm2notation);
    Assert.assertEquals(formula, "C18H32O16");
  }

  @Test
  public void testDirectNonIntegerRepeatIsRejectedAtBuildTime() throws ParserException, JDOMException,
      NotationException, ChemistryException, CTKException, PolymerIDsException, MonomerException,
      GroupingNotationException, ConnectionNotationException, MonomerLoadingException,
      org.helm.notation2.parser.exceptionparser.NotationException {
    String notation = "CARB1{[b-D-Glcp]'n'}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    try {
      MoleculePropertyCalculator.getMolecularFormular(helm2notation);
      Assert.fail("Expected BuilderMoleculeException for a non-integer direct CARB count");
    } catch (BuilderMoleculeException e) {
      Assert.assertTrue(e.getMessage().contains("non-integer count"), e.getMessage());
    }
  }

  @Test
  public void testNonIntegerRepeatIsRejectedAtBuildTime() throws ParserException, JDOMException,
      NotationException, ChemistryException, CTKException, PolymerIDsException, MonomerException,
      GroupingNotationException, ConnectionNotationException, MonomerLoadingException,
      org.helm.notation2.parser.exceptionparser.NotationException {
    // An open-ended repeat count ("n") has no concrete structure to build. Like the
    // other polymer types' non-integer repeats, it must parse and validate but fail
    // cleanly at build with a dedicated message, not silently mis-build.
    String notation = "CARB1{[b-D-Glcp].(R4:[b-D-Glcp])'n'}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    try {
      MoleculePropertyCalculator.getMolecularFormular(helm2notation);
      Assert.fail("Expected BuilderMoleculeException for a non-integer repeat count");
    } catch (BuilderMoleculeException e) {
      Assert.assertTrue(e.getMessage().contains("non-integer count"), e.getMessage());
    }
  }

  @Test
  public void testUnknownMonomerIsRejectedAtBuildTime() throws ParserException, JDOMException,
      NotationException, ChemistryException, CTKException, PolymerIDsException, MonomerException,
      GroupingNotationException, ConnectionNotationException, MonomerLoadingException,
      org.helm.notation2.parser.exceptionparser.NotationException {
    // A fully-unknown monomer has no concrete structure - it must validate but
    // fail cleanly at build, not silently mis-build.
    String notation = "CARB1{[b-D-Glcp].X}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    try {
      MoleculePropertyCalculator.getMolecularFormular(helm2notation);
      Assert.fail("Expected BuilderMoleculeException for an unknown CARB monomer");
    } catch (BuilderMoleculeException e) {
      Assert.assertTrue(e.getMessage().contains("Unknown CARB monomer"), e.getMessage());
    }
  }

  @Test
  public void testAmbiguityGroupIsRejectedAtBuildTime() throws ParserException, JDOMException,
      NotationException, ChemistryException, CTKException, PolymerIDsException, MonomerException,
      GroupingNotationException, ConnectionNotationException, MonomerLoadingException,
      org.helm.notation2.parser.exceptionparser.NotationException {
    // An ambiguity mixture cannot be resolved to a single molecule - build must
    // fail cleanly with a dedicated message.
    String notation = "CARB1{([b-D-Glcp]+[b-D-Glcp])}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
    try {
      MoleculePropertyCalculator.getMolecularFormular(helm2notation);
      Assert.fail("Expected BuilderMoleculeException for a CARB ambiguity group");
    } catch (BuilderMoleculeException e) {
      Assert.assertTrue(e.getMessage().contains("ambiguity group"), e.getMessage());
    }
  }

  /**
   * beta-D-glucopyranose, with R-group markers at every free hydroxyl per the
   * carbon numbering used in the notation (R1 anomeric; R2/R3/R4 ring
   * hydroxyls; R6 the exocyclic primary alcohol; no R5, since C5 bears no
   * free hydroxyl).
   */
  private Monomer glucopyranose() {
    Monomer m = new Monomer();
    m.setAlternateId("b-D-Glcp");
    m.setPolymerType(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    m.setMonomerType(Monomer.BACKBONE_MOMONER_TYPE);
    m.setName("beta-D-Glucopyranose");
    m.setCanSMILES("[H:6]OC[C@H]1O[C@@H]([OH:1])[C@H](O[H:2])[C@@H](O[H:3])[C@@H]1O[H:4]");
    m.getAttachmentList().add(attachment("R1", "OH", "O[*:1]"));
    m.getAttachmentList().add(attachment("R2", "H", "[*:2][H]"));
    m.getAttachmentList().add(attachment("R3", "H", "[*:3][H]"));
    m.getAttachmentList().add(attachment("R4", "H", "[*:4][H]"));
    m.getAttachmentList().add(attachment("R6", "H", "[*:6][H]"));
    return m;
  }

  private Attachment attachment(String label, String capGroupName, String capGroupSMILES) {
    Attachment att = new Attachment(label, capGroupName);
    att.setAlternateId(label + "-" + capGroupName);
    att.setCapGroupSMILES(capGroupSMILES);
    return att;
  }

}
