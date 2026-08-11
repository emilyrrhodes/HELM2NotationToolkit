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
    MonomerFactory.getInstance().getMonomerStore().addMonomer(glucopyranose());
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

  /**
   * beta-D-glucopyranose, with R-group markers at every free hydroxyl per the
   * carbon numbering used in the notation (R1 anomeric; R2/R3/R4 ring
   * hydroxyls; R6 the exocyclic primary alcohol; no R5, since C5 bears no
   * free hydroxyl).
   */
  private Monomer glucopyranose() {
    Monomer m = new Monomer();
    m.setAlternateId("Glcp");
    m.setPolymerType(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    m.setMonomerType(Monomer.BACKBONE_MOMONER_TYPE);
    m.setName("D-Glucopyranose");
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
