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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.helm.chemtoolkit.CTKException;
import org.helm.notation2.Chemistry;
import org.helm.notation2.InterConnections;
import org.helm.notation2.Monomer;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.MonomerStore;
import org.helm.notation2.exception.AttachmentException;
import org.helm.notation2.exception.ChemistryException;
import org.helm.notation2.exception.ConnectionNotationException;
import org.helm.notation2.exception.GroupingNotationException;
import org.helm.notation2.exception.HELM2HandledException;
import org.helm.notation2.exception.MonomerException;
import org.helm.notation2.exception.MonomerLoadingException;
import org.helm.notation2.exception.NotationException;
import org.helm.notation2.exception.PolymerIDsException;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.parser.notation.ValidationMethod;
import org.helm.notation2.parser.notation.connection.ConnectionNotation;
import org.helm.notation2.parser.notation.grouping.GroupingElement;
import org.helm.notation2.parser.notation.grouping.GroupingNotation;
import org.helm.notation2.parser.notation.polymer.GroupEntity;
import org.helm.notation2.parser.notation.polymer.HELMEntity;
import org.helm.notation2.parser.notation.polymer.MonomerNotation;
import org.helm.notation2.parser.notation.polymer.MonomerNotationGroup;
import org.helm.notation2.parser.notation.polymer.MonomerNotationGroupElement;
import org.helm.notation2.parser.notation.polymer.MonomerNotationList;
import org.helm.notation2.parser.notation.polymer.MonomerNotationUnit;
import org.helm.notation2.parser.notation.polymer.CarbBranch;
import org.helm.notation2.parser.notation.polymer.CarbMonomerNotationUnit;
import org.helm.notation2.parser.notation.polymer.CarbMonomerParser;
import org.helm.notation2.parser.notation.polymer.MonomerNotationUnitRNA;
import org.helm.notation2.parser.notation.polymer.PolymerNotation;
import org.helm.notation2.parser.notation.polymer.RNAEntity;
import org.jdom2.JDOMException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Validation class to validate the whole HELM string
 *
 * @author hecht
 */
public final class Validation {

	private static final Logger LOG = LoggerFactory.getLogger(Validation.class);

	/**
	 * Default constructor.
	 */
	private Validation() {

	}

	/**
	 * method to check if the generated notation objects by the parser are
	 * correct the polymer ids have to be unique; all monomers have to be valid;
	 * all used polymer ids in the grouping section have to be there; all
	 * connections have to be valid
	 *
	 * @param helm2notation
	 *            HELM2Notation object
	 * @throws PolymerIDsException
	 *             if the polymer section is not valid
	 * @throws MonomerException
	 *             if a monomer is not valid
	 * @throws GroupingNotationException
	 *             if the grouping section is not valid
	 * @throws ConnectionNotationException
	 *             if the connection section is not valid
	 * @throws NotationException
	 *             if notation is not valid
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 * @throws MonomerLoadingException
	 *             monomers can not be loaded
	 * @throws org.helm.notation2.parser.exceptionparser.NotationException
	 *             if notation is not valid
	 */
	public static void validateNotationObjects(HELM2Notation helm2notation) throws PolymerIDsException,
			MonomerException, GroupingNotationException, ConnectionNotationException, NotationException,
			ChemistryException, MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
		LOG.info("Validation process is starting");
		/* all polymer ids have to be unique */
		if (!validateUniquePolymerIDs(helm2notation)) {
			LOG.info("Polymer IDS have to be unique");
			throw new PolymerIDsException("Polymer IDs have to be unique");
		}
		/* Validation of Monomers */
		if (!validateMonomers(MethodsMonomerUtils.getListOfMonomerNotation(helm2notation.getListOfPolymers()))) {
			LOG.info("Monomers have to be valid");
			throw new MonomerException("Monomers have to be valid");
		}
		/* validate the grouping section */
		if (!validateGrouping(helm2notation)) {
			LOG.info("Group information is not valid");
			throw new GroupingNotationException("Group notation is not valid");
		}
		/* validate the connection */
		if (!validateConnections(helm2notation)) {
			LOG.info("Connection information is not valid");
			throw new ConnectionNotationException("Connection notation is not valid");
		}
	}

	/**
	 * method to validate a list of MonomerNotation objects
	 *
	 * @param mon
	 *            List of MonomerNotation objects
	 * @return true if all monomers are valid, false otherwise
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 * @throws MonomerLoadingException
	 *             if monomer loading was not successfull
	 * @throws org.helm.notation2.parser.exceptionparser.NotationException
	 *             if notation is not valid
	 */
	protected static boolean validateMonomers(List<MonomerNotation> mon) throws ChemistryException,
			MonomerLoadingException, org.helm.notation2.parser.exceptionparser.NotationException {
		MonomerStore monomerStore = MonomerFactory.getInstance().getMonomerStore();
		for (MonomerNotation monomerNotation : mon) {
			if (Monomer.CARBOHYDRATE_POLYMER_TYPE.equals(monomerNotation.getType())) {
				if (!validateCarbElement(monomerNotation, monomerStore)) {
					LOG.info("CARB monomer is not valid: " + monomerNotation.getUnit());
					return false;
				}
				continue;
			}
			if (!(isMonomerValid(monomerNotation.getUnit(), monomerNotation.getType()))) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Validates a single CARB polymer element. A concrete monomer must resolve to a
	 * base sugar in the store and reference only defined attachment points; a
	 * fully-unknown monomer ("X"/"*"/"?") is accepted verbatim (it carries no
	 * connection points); an ambiguity group (mixture "(A+B)" / or-group "(A,B)")
	 * is valid when every one of its member monomers is valid.
	 *
	 * @param monomerNotation the CARB element to validate
	 * @param monomerStore the monomer store to resolve base sugars against
	 * @return true if the element is valid
	 */
	private static boolean validateCarbElement(MonomerNotation monomerNotation, MonomerStore monomerStore) {
		if (monomerNotation instanceof CarbMonomerNotationUnit) {
			CarbMonomerNotationUnit carbUnit = (CarbMonomerNotationUnit) monomerNotation;
			/* an unknown monomer carries no connection points - nothing to resolve or check */
			if (carbUnit.isUnknown()) {
				return true;
			}
			Monomer monomer = resolveCarbMonomer(carbUnit.getUnit(), monomerStore);
			return monomer != null && isCarbAttachmentValid(carbUnit, monomer);
		}
		if (monomerNotation instanceof MonomerNotationGroup) {
			for (MonomerNotationGroupElement element : ((MonomerNotationGroup) monomerNotation).getListOfElements()) {
				if (!validateCarbElement(element.getMonomerNotation(), monomerStore)) {
					return false;
				}
			}
			return true;
		}
		/* a bare "?"/"_" unknown token still routes through the generic whitelist */
		return isMonomerValidUnknownToken(monomerNotation.getUnit());
	}

	/**
	 * @param str a monomer token
	 * @return true if the token is one of the generic unknown-monomer tokens ("?"/"_")
	 */
	private static boolean isMonomerValidUnknownToken(String str) {
		return str.equals("?") || str.equals("_") || str.equals("*") || str.equals("X");
	}

	/**
	 * method to check that every R-group a CARB monomer notation references -
	 * the R-group it is attached to its predecessor with, its own anomeric
	 * R-group, and the convergence R-group of each of its branches - is an
	 * attachment point that actually exists on the resolved monomer
	 *
	 * @param unit CARB monomer notation to check
	 * @param monomer the already-resolved base monomer
	 * @return true if every referenced R-group is a valid attachment point (or
	 *         unknown, "R?"), false otherwise
	 */
	private static boolean isCarbAttachmentValid(CarbMonomerNotationUnit unit, Monomer monomer) {
		if (!isCarbRGroupValid(unit.getIncomingRGroup(), monomer)
				|| !isCarbRGroupValid(unit.getAnomericRGroup(), monomer)) {
			return false;
		}
		for (CarbBranch branch : unit.getBranches()) {
			if (!isCarbRGroupValid(branch.getConvergenceRGroup(), monomer)) {
				return false;
			}
		}
		return true;
	}

	private static boolean isCarbRGroupValid(String rGroup, Monomer monomer) {
		if (rGroup == null || rGroup.equalsIgnoreCase("R?")) {
			return true;
		}
		return monomer.getAttachment(rGroup) != null;
	}

	/**
	 * method to resolve a CARB monomer to its monomer store entry
	 *
	 * <p>Monomers are keyed per stereoisomer: the full residue id as written -
	 * including the anomeric descriptor (a/b) and absolute configuration (D/L) -
	 * is the store key, so {@code a-D-Glcp} and {@code b-D-Glcp} resolve to
	 * distinct entries. The grammar guarantees canonical casing and the store's
	 * inner map is case-insensitive, so the raw id is a safe key. The residue is
	 * still parsed first, purely to reject malformed ids.
	 *
	 * @param carbId CARB monomer unit string, e.g. "[a-D-Glcp]" or "a-D-Glcp"
	 * @param monomerStore monomer store to resolve the monomer from
	 * @return the resolved monomer, or null if the unit does not parse or has
	 *         no matching entry in the store
	 */
	private static Monomer resolveCarbMonomer(String carbId, MonomerStore monomerStore) {
		if (carbId.startsWith("[") && carbId.endsWith("]")) {
			carbId = carbId.substring(1, carbId.length() - 1);
		}
		try {
			/* parse only to reject malformed residue content; the full id is the key */
			CarbMonomerParser.parse(carbId);
			return monomerStore.hasMonomer(Monomer.CARBOHYDRATE_POLYMER_TYPE, carbId)
					? monomerStore.getMonomer(Monomer.CARBOHYDRATE_POLYMER_TYPE, carbId) : null;
		} catch (org.helm.notation2.parser.exceptionparser.NotationException e) {
			return null;
		}
	}

	/**
	 * method to valid all existent connections in the Notation objects
	 *
	 * @param helm2notation
	 *            HELM2Notation object
	 * @return true if all connections are valid, false otherwise
	 * @throws NotationException
	 *             if notation is not valid
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 */
	public static boolean validateConnections(HELM2Notation helm2notation)
			throws NotationException, ChemistryException {
		try {
			LOG.info("Validation of Connection section starts");
			List<ConnectionNotation> listConnections = helm2notation.getListOfConnections();
			List<String> listPolymerIDs = helm2notation.getPolymerAndGroupingIDs();

			/* Hash-Map to save only specific InterConnections */
			InterConnections interconnection = new InterConnections();
			boolean specific = true;
			/* check for each single connection */
			for (ConnectionNotation connection : listConnections) {

				/* check polymer ids */
				checkPolymerIDSConnection(connection, listPolymerIDs);

				/* check for unspecific interaction */
				if (connection.getSourceId() instanceof GroupEntity
						|| connection.getTargetId() instanceof GroupEntity) {
					specific = false;
				}

				/* check Monomers:-> can be number */
				PolymerNotation source = helm2notation.getPolymerNotation(connection.getSourceId().getId());
				String sourceUnit = connection.getSourceUnit();
				PolymerNotation target = helm2notation.getPolymerNotation(connection.getTargetId().getId());
				String targetUnit = connection.getTargetUnit();

				/* check for specific interactions */
				if (isConnectionSpecific(connection) && specific) {
					/*
					 * interaction seems to be specific: it is given in number
					 * -> place of monomer
					 */
					/* Get Monomers */
					specific = true;
					int occurenceOne = Integer.parseInt(sourceUnit);
					int occurenceTwo = Integer.parseInt(targetUnit);

					/*
					 * if the monomers are a group or a list of monomers -> is
					 * it no more specific
					 */
					/* can the two form a connection */

					List<Monomer> listMonomersOne;

					listMonomersOne = getAllMonomers(source.getMonomerNotation(occurenceOne), occurenceOne);

					List<Monomer> listMonomersTwo = getAllMonomers(target.getMonomerNotation(occurenceTwo),
							occurenceTwo);

					/* check each single Attachment */
					checkAttachment(listMonomersOne, listMonomersTwo, connection, helm2notation, interconnection,
							specific);

				} /* Unspecific Interaction */ else {
					List<Integer> listMonomerOccurencesOne = getOccurencesOfMonomerNotation(sourceUnit,
							connection.getSourceId(), helm2notation);
					List<Integer> listMonomerOccurencesTwo = getOccurencesOfMonomerNotation(targetUnit,
							connection.getTargetId(), helm2notation);
					/* ? - section has to be included */
					if (listMonomerOccurencesOne.isEmpty()) {
						for (Integer occurenceTwo : listMonomerOccurencesTwo) {
							List<Monomer> listMonomersTwo = getAllMonomers(target.getMonomerNotation(occurenceTwo),
									occurenceTwo);
							checkSingleAttachment(listMonomersTwo, connection.getrGroupTarget(), helm2notation,
									connection, interconnection, connection.getTargetId().getId());
						}
					}
					for (Integer occurenceOne : listMonomerOccurencesOne) {
						/* get Monomers */
						List<Monomer> listMonomersOne = getAllMonomers(source.getMonomerNotation(occurenceOne),
								occurenceOne);
						checkSingleAttachment(listMonomersOne, connection.getrGroupSource(), helm2notation, connection,
								interconnection, connection.getSourceId().getId());
						/* check single attachment */
						for (Integer occurenceTwo : listMonomerOccurencesTwo) {
							List<Monomer> listMonomersTwo = getAllMonomers(target.getMonomerNotation(occurenceTwo),
									occurenceTwo);
							checkSingleAttachment(listMonomersTwo, connection.getrGroupTarget(), helm2notation,
									connection, interconnection, connection.getTargetId().getId());
							checkAttachment(listMonomersOne, listMonomersTwo, connection, helm2notation,
									interconnection, false);

						}
					}
				}
			}

			return true;
		} catch (PolymerIDsException | AttachmentException | HELM2HandledException | MonomerException | IOException
				| org.helm.notation2.parser.exceptionparser.NotationException | CTKException e) {
			e.printStackTrace();
			LOG.info(e.getMessage());
			return false;
		}
	}

	/**
	 * method to check if the given connection is specific
	 *
	 * @param connectionNotation
	 *            input ConnectionNotation
	 * @return true if the described connection is specific, false otherwise
	 */
	private static boolean isConnectionSpecific(ConnectionNotation connectionNotation) {
		String connection = connectionNotation.getSourceUnit() + ":" + connectionNotation.getrGroupSource() + "-"
				+ connectionNotation.getTargetUnit() + ":" + connectionNotation.getrGroupTarget();
		/* check for specific interaction */
		if (connection.matches("\\d+:R\\d-\\d+:R\\d|\\d+:pair-\\d+:pair")) {
			return true;
		}
		return false;
	}

	/**
	 * method to get all occurences of the MonomerNotation
	 *
	 * @param sourceUnit
	 *            given monomernotation to look for
	 * @param e
	 *            HELMEntity of the sourceUnit
	 * @param helm2notation
	 *            HELM2Notation object
	 * @return occurences of the MonomerNotation
	 * @throws org.helm.notation2.parser.exceptionparser.NotationException
	 * @throws AttachmentException
	 *             if attachment is not valid
	 */
	private static List<Integer> getOccurencesOfMonomerNotation(String sourceUnit, HELMEntity e,
			HELM2Notation helm2notation)
					throws org.helm.notation2.parser.exceptionparser.NotationException, AttachmentException {
		List<Integer> occurences = new ArrayList<Integer>();

		/* The monomer's position in the polymer is specified */
		try {
			occurences.add(Integer.parseInt(sourceUnit));
			return occurences;
		} catch (NumberFormatException ex) {
			MonomerNotation mon = ValidationMethod.decideWhichMonomerNotation(sourceUnit, e.getType());
			/* it is only one monomer e.g. C */
			if (mon instanceof MonomerNotationUnit) {
				PolymerNotation polymerNotation = helm2notation.getPolymerNotation(e.getId());
				/* monomer can also be unknown */
				if (sourceUnit.equals("?")) {
					return occurences;
				}

				occurences.addAll(findElementInPolymer(sourceUnit, polymerNotation));
				/* the specified monomer does not exist in the polymer */
				if (occurences.isEmpty()) {
					throw new AttachmentException("Monomer is not there");
				}
			} /* second: group (mixture or or) or list */ else if (mon instanceof MonomerNotationGroup
					|| mon instanceof MonomerNotationList) {
				PolymerNotation polymerNotation = helm2notation.getPolymerNotation(e.getId());
				Map<String, String> elements = new HashMap<String, String>();
				for (MonomerNotationGroupElement groupElement : ((MonomerNotationGroup) mon).getListOfElements()) {
					elements.put(groupElement.getMonomerNotation().getUnit(), "");
				}

				for (String e1 : elements.keySet()) {
					try {
						int i = Integer.parseInt(e1);
						elements.put(e1, "1");
						occurences.add(i);
					} catch (NumberFormatException ex1) {
						// have to be found in polymer
						List<Integer> foundMonomers = findElementInPolymer(e1, polymerNotation);
						if (foundMonomers.size() > 0) {
							elements.put(e1, "1");
							occurences.addAll(foundMonomers);
						}

					}
				}
				if (occurences.size() < elements.size() || elements.containsValue("")) {
					throw new AttachmentException("Not all Monomers are there");
				}
			}

			return occurences;
		}
	}

	private static List<Integer> findElementInPolymer(String e1, PolymerNotation polymerNotation) {
		List occurences = new ArrayList<Integer>();
		int j = 0;
		for (int i = 0; i < polymerNotation.getPolymerElements().getListOfElements().size(); i++) {
			MonomerNotation el = polymerNotation.getPolymerElements().getListOfElements().get(i);
			if (el instanceof MonomerNotationUnitRNA) {
				for (MonomerNotationUnit unit : ((MonomerNotationUnitRNA) el).getContents()) {
					j += 1;
					if (unit.getUnit().equals(e1)) {
						occurences.add(j);
					}

				}
			} else {
				j += 1;
				if (polymerNotation.getPolymerElements().getListOfElements().get(i).getUnit().equals(e1)) {
					occurences.add(i + 1);
				}
			}
		}

		return occurences;
	}

	/**
	 * method to validate every GroupNotation of the Notation objects
	 *
	 * @param helm2notation
	 *            HELM2Notation object
	 * @return true if the grouping is valid, false otherwise
	 */
	public static boolean validateGrouping(HELM2Notation helm2notation) {
		List<GroupingNotation> listGroupings = helm2notation.getListOfGroupings();
		List<String> listPolymerIDs = helm2notation.getPolymerAndGroupingIDs();

		/* validate each group */
		for (GroupingNotation grouping : listGroupings) {
			/* check for each group element if the polymer id is there */
			for (GroupingElement groupingElement : grouping.getAmbiguity().getListOfElements()) {
				if (!(listPolymerIDs.contains(groupingElement.getID().getId()))) {
					LOG.info("Element of Group: " + groupingElement.getID().getId() + " does not exist");
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * method to check if all existent polymer ids are unique
	 *
	 * @param helm2notation
	 *            HELM2Notation object
	 * @return true if all polymers are unique, false otherwise
	 */
	public static boolean validateUniquePolymerIDs(HELM2Notation helm2notation) {
		List<String> listPolymerIDs = helm2notation.getPolymerAndGroupingIDs();
		Map<String, String> uniqueId = new HashMap<String, String>();
		for (String polymerID : listPolymerIDs) {
			uniqueId.put(polymerID, "");
		}
		if (listPolymerIDs.size() > uniqueId.size()) {
			LOG.info("Polymer node IDs are not unique");
			return false;
		}

		return true;
	}

	/**
	 * method to check if the given polymer id exists in the given list of
	 * polymer ids
	 *
	 * @param str
	 *            polymer id
	 * @param listPolymerIDs
	 *            List of polymer ids
	 * @return true if the polymer id exists, false otherwise
	 * @throws PolymerIDsException
	 *             if the polymer id does not exist
	 */
	private static void checkExistenceOfPolymerID(String str, List<String> listPolymerIDs) throws PolymerIDsException {
		if (!(listPolymerIDs.contains(str))) {
			LOG.info("Polymer Id does not exist");
			throw new PolymerIDsException("Polymer ID does not exist");
		}
	}

	/**
	 * method to check the monomer's validation
	 *
	 * @param str
	 *            monomer id
	 * @param type
	 *            type of monomer
	 * @return true if the monomer is valid, false otherwise
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 * @throws MonomerLoadingException
	 *             if momomers can not be loaded
	 * @throws org.helm.notation2.parser.exceptionparser.NotationException
	 */
	private static boolean isMonomerValid(String str, String type) throws ChemistryException, MonomerLoadingException,
			org.helm.notation2.parser.exceptionparser.NotationException {
		LOG.info("Is Monomer valid: " + str);
		MonomerFactory monomerFactory = null;
		monomerFactory = MonomerFactory.getInstance();

		/* Search in Database */
		MonomerStore monomerStore = monomerFactory.getMonomerStore();
		if (monomerStore.hasMonomer(type, str)) {
			LOG.info("Monomer is located in the database: " + str);
			return true;
		} else if (str.charAt(0) == '[' && str.charAt(str.length() - 1) == ']'
				&& monomerStore.hasMonomer(type, str.substring(1, str.length() - 1))) {
			LOG.info("Monomer is located in the database: " + str);
			return true;
		} /* polymer type is Blob: accept all */ else if (type.equals("BLOB")) {
			LOG.info("Blob's Monomer Type: " + str);
			return true;
		} /* new unknown monomer for peptide */ else if (type.equals("PEPTIDE") && str.equals("X")) {
			LOG.info("Unknown monomer type for peptide: " + str);
			return true;
		} /* new unknown monomer for peptide */ else if (type.equals("RNA") && str.equals("N")) {
			LOG.info("Unknown monomer type for rna: " + str);
			return true;
		} /* new unknown types */ else if (str.equals("?") || str.equals("_")) {
			LOG.info("Unknown types: " + str);
			return true;
		} /* nucleotide */ else if (type.equals("RNA")) {
			List<String> elements = NucleotideParser.getMonomerIDListFromNucleotide(str);
			for (String element : elements) {
				if (!(monomerStore.hasMonomer(type, element))) {
					/* SMILES Check */
					if (element.startsWith("[") && element.endsWith("]")) {
						element = element.substring(1, element.length() - 1);
					}
					if (!Chemistry.getInstance().getManipulator().validateSMILES(element)) {
						return false;
					}
				}
			}
			LOG.info("Nucleotide type for RNA: " + str);
			return true;

		}

		LOG.info("SMILES Check");
		/* SMILES Check */
		if (str.charAt(0) == '[' && str.charAt(str.length() - 1) == ']') {
			str = str.substring(1, str.length() - 1);
		}

		return Chemistry.getInstance().getManipulator().validateSMILES(str);

	}

	/**
	 * method to check for one connection if the two polymer ids exist
	 *
	 * @param not
	 *            ConnectionNotation
	 * @param listPolymerIDs
	 *            List of polymer ids
	 * @throws PolymerIDsException
	 *             if the polymer ids do not exist
	 */
	private static void checkPolymerIDSConnection(ConnectionNotation not, List<String> listPolymerIDs)
			throws PolymerIDsException {
		/* the polymer ids have to be there */
		checkExistenceOfPolymerID(not.getSourceId().getId(), listPolymerIDs);
		checkExistenceOfPolymerID(not.getTargetId().getId(), listPolymerIDs);
	}

	/**
	 * method to get for one MonomerNotation all valid contained monomers
	 *
	 *
	 * @param not
	 *            MonomerNotation
	 * @param position
	 *            position of monomer notation
	 * @return List of Monomer
	 * @throws HELM2HandledException
	 *             if HELM2 features were there
	 * @throws MonomerException
	 *             if monomer is not valid
	 * @throws NotationException
	 *             if notation is not valid
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 * @throws CTKException
	 *             general ChemToolKit exception passed to HELMToolKit
	 * @throws MonomerLoadingException
	 *             if monomer store can not be loaded
	 */
	public static List<Monomer> getAllMonomers(MonomerNotation not, int position) throws HELM2HandledException,
			MonomerException, NotationException, ChemistryException, CTKException, MonomerLoadingException {
		List<Monomer> monomers = new ArrayList<Monomer>();

		MonomerFactory monomerFactory = MonomerFactory.getInstance();
		MonomerStore monomerStore = monomerFactory.getMonomerStore();
		if (not instanceof MonomerNotationUnitRNA) {
			monomers.addAll(getMonomersRNA((MonomerNotationUnitRNA) not, monomerStore, position));

		} else if (not instanceof CarbMonomerNotationUnit) {
			/*
			 * CARB monomers are keyed per stereoisomer by their full residue id
			 * (e.g. "a-D-Glcp"), so resolve via the same helper used for
			 * validation rather than MethodsMonomerUtils.getMonomer - the helper
			 * strips brackets and validates the residue before the store lookup.
			 */
			Monomer monomer = resolveCarbMonomer(not.getUnit(), monomerStore);
			if (monomer == null) {
				throw new MonomerException("Defined Monomer is not in the database: " + not.getUnit());
			}
			monomers.add(monomer);
		} else if (not instanceof MonomerNotationUnit) {
			String id = not.getUnit();
			// if (id.startsWith("[") && id.endsWith("]")) {
			// id = id.substring(1, id.length() - 1);
			// }
			monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
		} else if (not instanceof MonomerNotationGroup) {
			for (MonomerNotationGroupElement groupElement : ((MonomerNotationGroup) not).getListOfElements()) {
				String id = groupElement.getMonomerNotation().getUnit();
				/*
				 * if (id.startsWith("[") && id.endsWith("]")) { id =
				 * id.substring(1, id.length() - 1); }
				 */
				monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
			}
		} else if (not instanceof MonomerNotationList) {
			for (MonomerNotation listElement : ((MonomerNotationList) not).getListofMonomerUnits()) {
				if (listElement instanceof MonomerNotationUnitRNA) {
					monomers.addAll(getMonomersRNA(((MonomerNotationUnitRNA) listElement), monomerStore, position));
				} else {
					String id = listElement.getUnit();
					/*
					 * if (id.startsWith("[") && id.endsWith("]")) { id =
					 * id.substring(1, id.length() - 1); }
					 */
					monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
				}
			}

		}

		return monomers;
	}

	/**
	 * method to get for one MonomerNotation all valid contained monomers. But
	 * only the Base
	 *
	 *
	 * @param not
	 *            MonomerNotation
	 * @return List of all base monomers
	 * @throws HELM2HandledException
	 *             if HELM2 features were there
	 * @throws MonomerException
	 *             if monomer is not valid
	 * @throws NotationException
	 *             if notation is not valid
	 * @throws ChemistryException
	 *             if the Chemistry Engine can not be initialized
	 * @throws CTKException
	 *             general ChemToolKit exception passed to HELMToolKit
	 * @throws MonomerLoadingException
	 *             if monomers can not be loaded
	 */
	public static List<Monomer> getAllMonomersOnlyBase(MonomerNotation not) throws HELM2HandledException,
			MonomerException, NotationException, ChemistryException, CTKException, MonomerLoadingException {
		LOG.debug("Get base for " + not);
		List<Monomer> monomers = new ArrayList<Monomer>();

		MonomerFactory monomerFactory = MonomerFactory.getInstance();
		MonomerStore monomerStore = monomerFactory.getMonomerStore();
		LOG.debug("Which MonomerNotationType " + not.getClass());
		if (not instanceof MonomerNotationUnitRNA) {
			LOG.debug("MonomerNotationUnitRNA");
			monomers.addAll(getMonomersRNAOnlyBase((MonomerNotationUnitRNA) not, monomerStore));

		} else if (not instanceof MonomerNotationUnit) {
			String id = not.getUnit();
			if (id.startsWith("[") && id.endsWith("]")) {
				id = id.substring(1, id.length() - 1);
			}
			monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
		} else if (not instanceof MonomerNotationGroup) {
			LOG.debug("MonomerNotationGroup");
			for (MonomerNotationGroupElement groupElement : ((MonomerNotationGroup) not).getListOfElements()) {
				String id = groupElement.getMonomerNotation().getUnit();
				if (id.startsWith("[") && id.endsWith("]")) {
					id = id.substring(1, id.length() - 1);
				}
				monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
			}
		} else if (not instanceof MonomerNotationList) {
			LOG.debug("MonomerNotationList");
			for (MonomerNotation listElement : ((MonomerNotationList) not).getListofMonomerUnits()) {
				if (listElement instanceof MonomerNotationUnitRNA) {
					monomers.addAll(getMonomersRNAOnlyBase(((MonomerNotationUnitRNA) listElement), monomerStore));
				} else {
					String id = listElement.getUnit();
					if (id.startsWith("[") && id.endsWith("]")) {
						id = id.substring(1, id.length() - 1);
					}
					monomers.add(MethodsMonomerUtils.getMonomer(not.getType(), id, ""));
				}
			}

		}

		return monomers;
	}

	/**
	 * method to check the attachment point's existence
	 *
	 * @param mon
	 *            Monomer
	 * @param str
	 *            Attachment point
	 * @throws AttachmentException
	 *             if the Attachment point is not there
	 */
	private static void checkAttachmentPoint(Monomer mon, String str) throws AttachmentException {
		if (!(mon.getAttachmentListString().contains(str))) {
			if (!(str.equals("?") || str.equals("pair")) && !mon.getAlternateId().equals("?")) {
				LOG.info("Attachment point for source is not there");
				throw new AttachmentException("Attachment point for source is not there: " + str);
			}
		}
	}

	/**
	 * method to check the validation of the attachment
	 *
	 * @param listMonomersOne
	 *            List of Monomers of the source
	 * @param listMonomersTwo
	 *            List of Monomers of the target
	 * @param not
	 *            ConnectionNotation
	 * @param helm2notation
	 *            HELM2Notation object
	 * @param interconnection
	 *            InterConnections
	 * @param spec
	 *            specificity of the connection
	 * @return true if it valid, false otherwise
	 * @throws AttachmentException
	 */
	private static void checkAttachment(List<Monomer> listMonomersOne, List<Monomer> listMonomersTwo,
			ConnectionNotation not, HELM2Notation helm2notation, InterConnections interconnection, boolean spec)
					throws AttachmentException {
		boolean specific = spec;
		if (listMonomersOne.size() > 1 || listMonomersTwo.size() > 1) {
			specific = false;
		}
		for (Monomer monomerOne : listMonomersOne) {
			for (Monomer monomerTwo : listMonomersTwo) {
				
				if(monomerOne.getCanSMILES().equals("*|X|N") || monomerTwo.getCanSMILES().equals("*|X|N")){
				throw new AttachmentException("Monomer (canonicalsmiles = *) is undefined and should be defined")	;
				}
				/* Rna-Basepair-hydrogen bonds */
				if (monomerOne.getPolymerType().equals("RNA") && monomerTwo.getPolymerType().equals("RNA")
						&& not.getrGroupSource().equals("pair") && not.getrGroupTarget().equals("pair")) {
					LOG.info("RNA strand connection");

					if (!(monomerOne.getMonomerType().equals("Branch")
							| monomerTwo.getMonomerType().equals("Branch"))) {
						LOG.info("RNA strand connection is not valid");
						throw new AttachmentException("RNA strand connection is not valid");
					}
			

					/*
					 * is the attachment point already occupied by another
					 * monomer?
					 */
					String detailsource = not.getSourceUnit() + "$" + not.getrGroupSource();
					String detailtarget = not.getTargetUnit() + "$" + not.getrGroupTarget();

					/*
					 * Is the attachment point already occupied by another
					 * monomer
					 */
					/* Intra connections */
					if (helm2notation.getSimplePolymer(not.getSourceId().getId()).getMapIntraConnection()
							.containsKey(detailsource)) {
						LOG.info("Attachment point is already occupied");
						throw new AttachmentException("Attachment point is already occupied");
					}
					if (helm2notation.getSimplePolymer(not.getTargetId().getId()).getMapIntraConnection()
							.containsKey(detailtarget)) {
						LOG.info("Attachment point is already occupied");
						throw new AttachmentException("Attachment point is already occupied");
					}
				}

				String detailsource = not.getSourceUnit() + "$" + not.getrGroupSource();
				String detailtarget = not.getTargetUnit() + "$" + not.getrGroupTarget();

				/* Inter connections */
				detailsource = not.getSourceId().getId() + "$" + detailsource;

				detailtarget = not.getTargetId().getId() + "$" + detailtarget;

				/* check */
				if (interconnection.hasKey(detailsource)) {
					LOG.info("Attachment point is already occupied");
					throw new AttachmentException("Attachment point is already occupied");
				}

				if (interconnection.hasKey(detailtarget)) {
					LOG.info("Attachment point is already occupied");
					throw new AttachmentException("Attachment point is already occupied");
				}
				/* Inter connections */
				detailsource = not.getSourceId().getId() + "$" + not.getSourceUnit() + "$" + not.getrGroupSource();
				detailtarget = not.getTargetId().getId() + "$" + not.getTargetUnit() + "$" + not.getrGroupTarget();

				if (specific) {
					/* save only specific interactions */
					interconnection.addConnection(detailsource, "");
					interconnection.addConnection(detailtarget, "");
				}
			}
		}
	}

	/**
	 * method to check for one attachment point the validation
	 *
	 * @param monomers
	 *            List of monomers
	 * @param rGroup
	 *            rGroup of the connection
	 * @param helm2notation
	 *            HELM2Notation object
	 * @param not
	 *            ConnectionNotation
	 * @param interconnection
	 *            InterConnections
	 * @param id
	 * @return true if it is valid, false otherwise
	 * @throws AttachmentException
	 */
	private static boolean checkSingleAttachment(List<Monomer> monomers, String rGroup, HELM2Notation helm2notation,
			ConnectionNotation not, InterConnections interconnection, String id) throws AttachmentException {

		for (Monomer monomer : monomers) {
			/* Are the attachment points there */
			checkAttachmentPoint(monomer, rGroup);

			/* is the attachment point already occupied by another monomer? */
			String detail = not.getSourceUnit() + "$" + not.getrGroupSource();

			/* Is the attachment point already occupied by another monomer */
			/* Intra connections */
			if (helm2notation.getSimplePolymer(id).getMapIntraConnection().containsKey(detail)) {
				throw new AttachmentException("Attachment point is already occupied");
			}

			/* Inter connections */
			detail = id + "$" + detail;

			/* check */
			if (interconnection.hasKey(detail)) {
				throw new AttachmentException("Attachment point is already occupied");
			}

		}
		return true;
	}

	/**
	 * method to get all monomers for MonomerNotationUnitRNA
	 *
	 * @param rna
	 *            MonomerNotationUnitRNA
	 * @param monomerStore
	 *            MonomerStore
	 * @return List of monomers of the MonomerNotationUnitRNA
	 * @throws HELM2HandledException
	 *             if HELM2 features were there
	 */
	private static List<Monomer> getMonomersRNA(MonomerNotationUnitRNA rna, MonomerStore monomerStore, int position)
			throws HELM2HandledException {
		try {
			List<Monomer> monomers = new ArrayList<Monomer>();
			for (int index = 0; index < rna.getContents().size(); index++) {
				String id = rna.getContents().get(index).getUnit();
				if (rna.getContents().get(index).getUnit().startsWith("[")
						&& rna.getContents().get(index).getUnit().endsWith("]")) {
					id = id.substring(1, id.length() - 1);
				}
				/* Special case */
				if (rna.getContents().size() == 1 && position == 0) {
					monomers.add(MethodsMonomerUtils.getMonomer(rna.getType(), id, "P"));

				} else {

					monomers.add(MethodsMonomerUtils.getMonomer(rna.getType(), id, rna.getInformation().get(index)));
				}
			}
			return monomers;
		} catch (Exception e) {
			e.printStackTrace();
			throw new HELM2HandledException(e.getMessage());
		}

	}

	/**
	 * method to get only the nucleotide base for one MonomerNotationUnitRNA
	 *
	 * @param rna
	 *            MonomerNotationUnitRNA
	 * @param monomerStore
	 *            MonomerStore
	 * @return List of Monomers
	 * @throws HELM2HandledException
	 *             if HELM2 features were there
	 */
	private static List<Monomer> getMonomersRNAOnlyBase(MonomerNotationUnitRNA rna, MonomerStore monomerStore)
			throws HELM2HandledException {
		try {
			List<Monomer> monomers = new ArrayList<Monomer>();
			for (MonomerNotationUnit unit : rna.getContents()) {
				String id = unit.getUnit().replace("[", "");
				id = id.replace("]", "");
				Monomer mon = MethodsMonomerUtils.getMonomer(rna.getType(), id, "");

				if (mon.getMonomerType().equals(Monomer.BRANCH_MOMONER_TYPE)) {
					monomers.add(mon);
				}
			}
			return monomers;
		} catch (Exception e) {
			e.printStackTrace();
			throw new HELM2HandledException(e.getMessage());
		}

	}
}
