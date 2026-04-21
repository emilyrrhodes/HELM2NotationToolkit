/**
 * *****************************************************************************
 * Copyright C 2015, The Pistoia Alliance
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *****************************************************************************
 */
package org.helm.notation2.tools;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.helm.chemtoolkit.CTKException;
import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.exception.ChemistryException;
import org.helm.notation2.exception.ConnectionNotationException;
import org.helm.notation2.exception.GroupingNotationException;
import org.helm.notation2.exception.HELM2HandledException;
import org.helm.notation2.exception.MonomerException;
import org.helm.notation2.exception.MonomerLoadingException;
import org.helm.notation2.exception.NotationException;
import org.helm.notation2.exception.ParserException;
import org.helm.notation2.exception.PolymerIDsException;
import org.helm.notation2.exception.RNAUtilsException;
import org.helm.notation2.parser.exceptionparser.HELM1ConverterException;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.tools.ChangeObjects;
import org.helm.notation2.tools.HELM2NotationUtils;
import org.jdom2.JDOMException;
import org.helm.notation2.parser.notation.polymer.CarbEntity;
import org.testng.Assert;
import org.testng.annotations.Test;

public class HELM2NotationUtilsTest {

  @Test
  public void testCARBPolymer() throws ParserException, JDOMException {
    String notation = "CARB1{[β-D-Gal].R4[β-D-GlcNAc].R6[α-D-GalNAc]}$$$$V2.0";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    Assert.assertEquals(helm2notation.getListOfPolymers().size(), 1);
    Assert.assertTrue(helm2notation.getListOfPolymers().get(0).getPolymerID() instanceof CarbEntity);
    Assert.assertEquals(helm2notation.getListOfPolymers().get(0).getPolymerElements().getListOfElements().size(), 3);
  }


	@Test
	public void testaddNotation()
			throws ParserException, JDOMException, org.helm.notation2.parser.exceptionparser.NotationException {
		int number = 0;
		HELM2Notation current = HELM2NotationUtils.readNotation(
				"PEPTIDE1{*}\"LC\"|PEPTIDE2{*}\"HC\"|PEPTIDE3{*}\"HC\"|PEPTIDE4{*}\"LC\"|CHEM1{*}$G1,CHEM1,K:R3-1:R1|PEPTIDE2,PEPTIDE3,250:R3-250:R3\"Hinge S-S connection\"|PEPTIDE2,PEPTIDE3,252:R3-252:R3\"Hinge S-Sconnection\"|PEPTIDE1,PEPTIDE2,120:R3-248:R3\"LC Hinge S-Sconnection\"|PEPTIDE4,PEPTIDE3,120:R3-248:R3\"LC Hinge S-Sconnection\"$G1(PEPTIDE1+PEPTIDE2+PEPTIDE3+PEPTIDE4)|G2(G1+CHEM1:4.5)$Hallo$V2.0");

		HELM2Notation toadd = HELM2NotationUtils.readNotation(
				"PEPTIDE1{A.X.G.C.(_,N).(A:10,G:30,R:30).T.C.F.D.W\"mutation\".(A:?+G:1.5).C}|RNA1{R(A)P.(R(N)P)'4'.(R(G)P)'3-7'\"mutation\"}|CHEM1{*}|BLOB1{BEAD}\"Animated Polystyrene\"$PEPTIDE1,BLOB1,X:R3-?:?\"Specific Conjugation\"|PEPTIDE1,CHEM1,(A+T):R3-?:?|PEPTIDE1,PEPTIDE1,(4,8):pair-12:pair$G1(PEPTIDE1:1+RNA1:2.5-2.7+BLOB1)|G2(G1:45,CHEM1:55)${\"Name\":\"lipid nanoparticle with RNA payload and peptide ligand\"}$V2.0");
		number = toadd.getListOfPolymers().size() + current.getListOfPolymers().size();
		HELM2NotationUtils.combineHELM2notation(current, toadd);

		Assert.assertEquals(current.getListOfPolymers().size(), number);
	}

	@Test
	public void testCount() throws ParserException, JDOMException, NotationException {
		HELM2Notation helm2notation = HELM2NotationUtils.readNotation("RNA1{R(A)P.[mR](A)}$$$$V2.0");
		Assert.assertEquals(HELM2NotationUtils.getTotalMonomerCount(helm2notation), 5);
	}

	@Test
	public void testReplaceMonomer() throws NotationException, MonomerLoadingException, MonomerException, JDOMException,
			IOException, ParserException, org.helm.notation2.parser.exceptionparser.NotationException,
			HELM1ConverterException, RNAUtilsException, HELM2HandledException, ChemistryException, CTKException {

		String notation = "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$";

		HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
		ChangeObjects.replaceMonomer(helm2notation, "RNA", "P", "sP");

		Assert.assertEquals(helm2notation.toHELM2(),
				"RNA1{R(A)[sP].R(G)[sP].R(C)[sP].R(U)[sP].R(A)[sP].R(A)[sP].R(A)[sP].R(G)[sP].R(G)}|RNA2{R(C)[sP].R(C)[sP].R(U)[sP].R(U)[sP].R(U)[sP].R(A)[sP].R(G)[sP].R(C)[sP].R(U)}$$$$V2.0");

		/**/
		notation = "PEPTIDE1{A'23'.C.D'12'.E'24'}|PEPTIDE2{G'22'.C.S'8'.P.P.P.P.P.P.P.P.P.K'6'}$$$$V2.0";

		helm2notation = HELM2NotationUtils.readNotation(notation);
		ChangeObjects.replaceMonomer(helm2notation, "PEPTIDE", "A", "G");

		Assert.assertEquals(helm2notation.toHELM2(),
				"PEPTIDE1{G'23'.C.D'12'.E'24'}|PEPTIDE2{G'22'.C.S'8'.P.P.P.P.P.P.P.P.P.K'6'}$$$$V2.0");

		notation = "RNA1{P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}|RNA2{P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$";
		helm2notation = HELM2NotationUtils.readNotation(notation);
		ChangeObjects.hybridize(helm2notation);
		System.out.println(helm2notation.toHELM2());
		Assert.assertEquals(helm2notation.getListOfConnections().size(), 6);
		notation = "RNA1{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}$$RNA2,RNA1,20:pair-2:pair|RNA2,RNA1,5:pair-17:pair|RNA2,RNA1,2:pair-20:pair|RNA2,RNA1,11:pair-11:pair|RNA2,RNA1,17:pair-5:pair|RNA2,RNA1,14:pair-8:pair|RNA2,RNA1,8:pair-14:pair$RNA1{ss}$";

		helm2notation = HELM2NotationUtils.readNotation(notation);
		String[] formatedSeqs = HELM2NotationUtils.getFormatedSirnaSequences(helm2notation);
		for (int i = 0; i < formatedSeqs.length; i++) {
			System.out.println(formatedSeqs[i]);
		}
	}

	@Test
  public void replaceSMILESIntoAlternatdeID() throws ParserException, JDOMException, PolymerIDsException, MonomerException, GroupingNotationException,
      ConnectionNotationException, NotationException, ChemistryException, org.helm.notation2.parser.exceptionparser.NotationException, IOException, HELM2HandledException, CTKException {
  

    String notation =
        "PEPTIDE1{([C[C@H](N[OH:1])C([H:2])=O]+G:?).G.G.G.C.C.K.K.K.K}|CHEM1{[MCC]}|RNA1{R(C)P.R([Nc1ncnc2n([*:1])cnc12])[sP].RP.R(G)P.[LR]([5meC])P}$PEPTIDE1,CHEM1,10:R3-1:R1$$$";
    HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
    ChangeObjects.replaceSMILESWithTemporaryIds(helm2notation);
    Pattern myPattern = Pattern.compile(("PM#[1-9][0-9]*"));

   String newNotation = helm2notation.toHELM2();
   Matcher matcher = myPattern.matcher(newNotation);
   String newMonomerID = "";
   while(matcher.find()){
	   newMonomerID = matcher.group();
   }
   
    
   Monomer m =  MonomerFactory.getInstance().getMonomerStore().getMonomer("PEPTIDE",newMonomerID);
   
   Assert.assertEquals(m.getAttachment("R1").getAlternateId(),"R1-OH");
   Assert.assertEquals(m.getAttachment("R2").getAlternateId(),"R2-H");

  }

}
