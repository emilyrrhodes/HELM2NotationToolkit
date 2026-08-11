package org.helm.notation2.tools;

import java.io.IOException;

import java.lang.reflect.Field;

import org.helm.notation2.Attachment;
import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.exception.ChemistryException;
import org.helm.notation2.wsadapter.MonomerStoreConfiguration;
import org.helm.notation2.exception.ConnectionNotationException;
import org.helm.notation2.exception.GroupingNotationException;
import org.helm.notation2.exception.HELM1FormatException;
import org.helm.notation2.exception.MonomerException;
import org.helm.notation2.exception.MonomerLoadingException;
import org.helm.notation2.exception.NotationException;
import org.helm.notation2.exception.ParserException;
import org.helm.notation2.exception.PolymerIDsException;
import org.helm.notation2.exception.ValidationException;
import org.helm.chemtoolkit.CTKException;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.parser.notation.polymer.CarbEntity;
import org.helm.notation2.parser.notation.polymer.CarbMonomerNotation;
import org.helm.notation2.parser.notation.polymer.CarbMonomerParser;
import org.jdom2.JDOMException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class CarbValidationTest {

  @BeforeClass
  public void registerCarbMonomers() throws IOException, MonomerException, ChemistryException,
      NoSuchFieldException, IllegalAccessException {
    Field wsField = MonomerStoreConfiguration.class.getDeclaredField("isUseWebservice");
    wsField.setAccessible(true);
    wsField.set(MonomerStoreConfiguration.getInstance(), false);
    // Monomers are keyed by base name; a/b and D/L are qualifiers parsed at lookup time.
    // addMonomer() silently skips a base name that is already present, so evict any
    // pre-existing (possibly stale, cache-persisted) CARB entry first to keep this
    // test hermetic regardless of the shared ~/.helm monomer cache state.
    registerCarb(carbMonomer("Gal", "Galactose"));
    registerCarb(carbMonomer("GlcNAc", "N-Acetylglucosamine"));
    registerCarb(carbMonomer("GalNAc", "N-Acetylgalactosamine"));
  }

  @Test
  public void testCarbMonomerParsing() throws org.helm.notation2.parser.exceptionparser.NotationException {
    CarbMonomerNotation gal = CarbMonomerParser.parse("b-D-Gal");
    Assert.assertEquals(gal.getAnomer(), "b");
    Assert.assertEquals(gal.getAbsoluteConfiguration(), "D");
    Assert.assertEquals(gal.getBaseName(), "Gal");

    CarbMonomerNotation galNAc = CarbMonomerParser.parse("a-D-GalNAc");
    Assert.assertEquals(galNAc.getAnomer(), "a");
    Assert.assertEquals(galNAc.getAbsoluteConfiguration(), "D");
    Assert.assertEquals(galNAc.getBaseName(), "GalNAc");

    CarbMonomerNotation unqualified = CarbMonomerParser.parse("Gal");
    Assert.assertNull(unqualified.getAnomer());
    Assert.assertNull(unqualified.getAbsoluteConfiguration());
    Assert.assertEquals(unqualified.getBaseName(), "Gal");
  }

  @Test
  public void testCARBPolymerIsRecognized() throws ParserException, JDOMException {
    String notation = "CARB1{[b-D-Gal].R4:[b-D-GlcNAc].R6:[a-D-GalNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Assert.assertEquals(helm2notation.getListOfPolymers().size(), 1);
    Assert.assertTrue(helm2notation.getListOfPolymers().get(0).getPolymerID() instanceof CarbEntity);
    Assert.assertEquals(helm2notation.getListOfPolymers().get(0).getPolymerElements().getListOfElements().size(), 3);
  }

  @Test
  public void testCARBValidation() throws ParserException, JDOMException,
      PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // b-D-Gal b(1->4) b-D-GlcNAc forms LacNAc; a-D-GalNAc initiates O-glycan core structures
    String notation = "CARB1{[b-D-Gal].R4:[b-D-GlcNAc].R6:[a-D-GalNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
  }

  @Test(expectedExceptions = MonomerException.class)
  public void testCARBValidationRejectsUndefinedRGroup() throws ParserException, JDOMException,
      PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // GlcNAc has no R9 attachment point defined - must be rejected, not silently accepted
    String notation = "CARB1{[b-D-Gal].R9:[b-D-GlcNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
  }

  @Test
  public void testCARBUnknownMonomerValidates() throws ParserException, JDOMException,
      PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // A fully-unknown monomer carries no connection points and needs no store
    // lookup - a lone unknown polymer, and an unknown mixed into a real chain,
    // must both pass validation (they only fail later, at build time).
    Validation.validateNotationObjects(HELM2NotationUtils.readNotation("CARB1{*}$$$$V2.0"));
    Validation.validateNotationObjects(
        HELM2NotationUtils.readNotation("CARB1{[b-D-Gal].X.R4:[b-D-GlcNAc]}$$$$V2.0"));
  }

  @Test
  public void testCARBAmbiguityGroupsValidate() throws ParserException, JDOMException,
      PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // A mixture and an or-group are valid when every member monomer resolves and
    // references only defined attachment points.
    Validation.validateNotationObjects(
        HELM2NotationUtils.readNotation("CARB1{([b-D-Gal]+[b-D-GlcNAc])}$$$$V2.0"));
    Validation.validateNotationObjects(
        HELM2NotationUtils.readNotation("CARB1{([b-D-Gal],[b-D-GlcNAc])}$$$$V2.0"));
  }

  @Test(expectedExceptions = MonomerException.class)
  public void testCARBAmbiguityGroupWithUndefinedMemberIsRejected() throws ParserException, JDOMException,
      PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException,
      MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
    // An or-group is only valid if every member is - here "Xyz" is not in the store.
    String notation = "CARB1{([b-D-Gal],[b-D-Xyz])}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
  }

  @Test(expectedExceptions = HELM1FormatException.class)
  public void testCARBIsRejectedByHELM1Conversion() throws ParserException, JDOMException, MonomerLoadingException,
      CTKException, ValidationException, ChemistryException, HELM1FormatException {
    // HELM1 predates glycans and has no CARB equivalent - conversion must fail
    // clearly rather than silently mishandle or drop the CARB polymer
    String notation = "CARB1{[b-D-Gal].R4:[b-D-GlcNAc].R6:[a-D-GalNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    HELM1Utils.getStandard(helm2notation);
  }

  private void registerCarb(Monomer monomer) throws IOException, MonomerException, ChemistryException {
    java.util.Map<String, java.util.Map<String, Monomer>> db =
        MonomerFactory.getInstance().getMonomerStore().getMonomerDB();
    java.util.Map<String, Monomer> carbMap = db.get(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    if (carbMap != null) {
      carbMap.remove(monomer.getAlternateId());
    }
    MonomerFactory.getInstance().getMonomerStore().addMonomer(monomer);
  }

  private Monomer carbMonomer(String alternateId, String name) {
    Monomer m = new Monomer();
    m.setAlternateId(alternateId);
    m.setPolymerType(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    m.setMonomerType(Monomer.BACKBONE_MOMONER_TYPE);
    m.setName(name);
    m.getAttachmentList().add(new Attachment("R1", Attachment.CAP_GROUP_H));
    m.getAttachmentList().add(new Attachment("R3", Attachment.CAP_GROUP_OH));
    m.getAttachmentList().add(new Attachment("R4", Attachment.CAP_GROUP_OH));
    m.getAttachmentList().add(new Attachment("R6", Attachment.CAP_GROUP_OH));
    return m;
  }

}
