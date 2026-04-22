package org.helm.notation2.tools;

import java.io.IOException;

import java.lang.reflect.Field;

import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.exception.ChemistryException;
import org.helm.notation2.wsadapter.MonomerStoreConfiguration;
import org.helm.notation2.exception.ConnectionNotationException;
import org.helm.notation2.exception.GroupingNotationException;
import org.helm.notation2.exception.MonomerException;
import org.helm.notation2.exception.MonomerLoadingException;
import org.helm.notation2.exception.NotationException;
import org.helm.notation2.exception.ParserException;
import org.helm.notation2.exception.PolymerIDsException;
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
    // Monomers are keyed by base name; a/b and D/L are qualifiers parsed at lookup time
    MonomerFactory.getInstance().getMonomerStore().addMonomer(carbMonomer("Gal", "Galactose"));
    MonomerFactory.getInstance().getMonomerStore().addMonomer(carbMonomer("GlcNAc", "N-Acetylglucosamine"));
    MonomerFactory.getInstance().getMonomerStore().addMonomer(carbMonomer("GalNAc", "N-Acetylgalactosamine"));
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
    String notation = "CARB1{[b-D-Gal].R4[b-D-GlcNAc].R6[a-D-GalNAc]}$$$$V2.0";
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
    // b-D-Gal b(1→4) b-D-GlcNAc forms LacNAc; a-D-GalNAc initiates O-glycan core structures
    String notation = "CARB1{[b-D-Gal].R4[b-D-GlcNAc].R6[a-D-GalNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Validation.validateNotationObjects(helm2notation);
  }

  private Monomer carbMonomer(String alternateId, String name) {
    Monomer m = new Monomer();
    m.setAlternateId(alternateId);
    m.setPolymerType(Monomer.CARBOHYDRATE_POLYMER_TYPE);
    m.setMonomerType(Monomer.BACKBONE_MOMONER_TYPE);
    m.setName(name);
    return m;
  }

}
