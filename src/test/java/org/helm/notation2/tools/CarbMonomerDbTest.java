package org.helm.notation2.tools;

import java.lang.reflect.Field;

import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.calculation.MoleculePropertyCalculator;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.wsadapter.MonomerStoreConfiguration;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Verifies the CARB monosaccharide monomers that ship in
 * {@code MonomerDBGZEncoded.xml}. Unlike {@link CarbBuilderMoleculeTest} and
 * {@link CarbValidationTest} (which register sugars programmatically), these
 * tests resolve monomers straight from the loaded, shipped database - proving
 * the new {@code <Polymer polymerType="CARB">} block loads, validates, and
 * builds chemically correct structures.
 */
public class CarbMonomerDbTest {

  @BeforeClass
  public void useLocalStore() throws NoSuchFieldException, IllegalAccessException {
    // Load monomers from the local/shipped store, not the web service (whose load
    // path does not handle CARB).
    Field wsField = MonomerStoreConfiguration.class.getDeclaredField("isUseWebservice");
    wsField.setAccessible(true);
    wsField.set(MonomerStoreConfiguration.getInstance(), false);
  }

  private String formula(String notation) throws Exception {
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    return MoleculePropertyCalculator.getMolecularFormular(helm2notation);
  }

  @Test
  public void testShippedMonomersBuildToKnownFormula() throws Exception {
    // A representative spread across the shape categories - hexose, N-acetyl-
    // hexosamine, deoxyhexose, pentose, uronic acid, and the sialic acid ketose -
    // each resolved from the shipped DB and built as a single free monomer.
    Assert.assertEquals(formula("CARB1{[b-D-Glcp]}$$$$V2.0"), "C6H12O6");
    Assert.assertEquals(formula("CARB1{[a-D-Manp]}$$$$V2.0"), "C6H12O6");
    Assert.assertEquals(formula("CARB1{[b-D-GlcpNAc]}$$$$V2.0"), "C8H15NO6");
    Assert.assertEquals(formula("CARB1{[a-L-Fucp]}$$$$V2.0"), "C6H12O5");
    Assert.assertEquals(formula("CARB1{[b-D-Xylp]}$$$$V2.0"), "C5H10O5");
    Assert.assertEquals(formula("CARB1{[b-D-GlcpA]}$$$$V2.0"), "C6H10O7");
    Assert.assertEquals(formula("CARB1{[a-D-Neup5Ac]}$$$$V2.0"), "C11H19NO9");
  }

  @Test
  public void testLacNAcDisaccharideLosesOneWater() throws Exception {
    // N-acetyllactosamine: b-D-Galp b(1->4) b-D-GlcpNAc. Gal (C6H12O6) + GlcNAc
    // (C8H15NO6) joined by one glycosidic bond, minus one H2O = C14H25NO11.
    Assert.assertEquals(formula("CARB1{[b-D-Galp].R4:[b-D-GlcpNAc]}$$$$V2.0"), "C14H25NO11");
  }

  @Test
  public void testAnomersAreDistinctEntries() throws Exception {
    // The per-stereoisomer keying means a and b anomers are separate monomers
    // with distinct (isomeric) SMILES - not one entry shared by base name.
    Monomer alpha = MonomerFactory.getInstance().getMonomerStore()
        .getMonomer(Monomer.CARBOHYDRATE_POLYMER_TYPE, "a-D-Glcp");
    Monomer beta = MonomerFactory.getInstance().getMonomerStore()
        .getMonomer(Monomer.CARBOHYDRATE_POLYMER_TYPE, "b-D-Glcp");
    Assert.assertNotNull(alpha, "a-D-Glcp should be in the shipped store");
    Assert.assertNotNull(beta, "b-D-Glcp should be in the shipped store");
    Assert.assertNotEquals(alpha.getCanSMILES(), beta.getCanSMILES(),
        "anomers must have distinct SMILES");
  }
}
