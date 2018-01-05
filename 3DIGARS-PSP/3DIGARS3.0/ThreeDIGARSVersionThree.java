import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ThreeDIGARSVersionThree {
	


	// alpha and beta values obtain from 3DIGARS
	// alpha-HP = 1.3802541
	// alpha-HH = 1.6832844
	// alpha-PP = 1.9315737
	// beta-HP = 1.4921875
	// beta-HH = 0.55859375
	// beta-PP = 0.265625
	// static float[] alphas = { 1.3802541f, 1.6832844f, 1.9315737f };
	static float[] betas = { 1.4921875f, 0.55859375f, 0.265625f };	// weight of contribution of HP, HH and PP
//	static float w1 = 0.5278592f; // weight for 3DIGARS2.0
	static float w1 = 0.7871094f, w2 = 0.009765625f, w3 = 0.013671875f;	// weights for 3DIGARS3.0

	

	// setup needed for loading frequency table and calculating energy
	static int maxRow = 14028;
	static int maxCol = 30;
	static double small_num = 0.000001;
	static int atomType = 167;

	// frequency table setup

	static float[][] alphaCarbonArray = new float[70000][3];
	static int[] resAtomPairArrayID = new int[70000];
	static String[] pairs = new String[maxRow];
	static int cAlphaCount = -1;
	static int[] resNum = new int[70000];
	static int[] Hphob_Hphil_Array = new int[70000];

	// frequency table setup complete

	static String[][] Hphob_Hphil_Mapper = { { "ALA", "hphob" },
			{ "ARG", "hphil" }, { "ASN", "hphil" }, { "ASP", "hphil" },
			{ "CYS", "hphil" }, { "GLN", "hphil" }, { "GLU", "hphil" },
			{ "GLY", "hphob" }, { "HIS", "hphil" }, { "ILE", "hphob" },
			{ "LEU", "hphob" }, { "LYS", "hphil" }, { "MET", "hphob" },
			{ "PHE", "hphob" }, { "PRO", "hphil" }, { "SER", "hphil" },
			{ "THR", "hphil" }, { "TRP", "hphil" }, { "TYR", "hphob" },
			{ "VAL", "hphob" } }; // this categorization is based on Dr. Hoque's
									// paper

	static String resAtomPair[] = { "ALA N", "ALA CA", "ALA C", "ALA O",
			"ALA CB", "CYS N", "CYS CA", "CYS C", "CYS O", "CYS CB", "CYS SG",
			"ASP N", "ASP CA", "ASP C", "ASP O", "ASP CB", "ASP CG", "ASP OD1",
			"ASP OD2", "GLU N", "GLU CA", "GLU C", "GLU O", "GLU CB", "GLU CG",
			"GLU CD", "GLU OE1", "GLU OE2", "PHE N", "PHE CA", "PHE C",
			"PHE O", "PHE CB", "PHE CG", "PHE CD1", "PHE CD2", "PHE CE1",
			"PHE CE2", "PHE CZ", "GLY N", "GLY CA", "GLY C", "GLY O", "HIS N",
			"HIS CA", "HIS C", "HIS O", "HIS CB", "HIS CG", "HIS ND1",
			"HIS CD2", "HIS CE1", "HIS NE2", "ILE N", "ILE CA", "ILE C",
			"ILE O", "ILE CB", "ILE CG1", "ILE CG2", "ILE CD1", "LYS N",
			"LYS CA", "LYS C", "LYS O", "LYS CB", "LYS CG", "LYS CD", "LYS CE",
			"LYS NZ", "LEU N", "LEU CA", "LEU C", "LEU O", "LEU CB", "LEU CG",
			"LEU CD1", "LEU CD2", "MET N", "MET CA", "MET C", "MET O",
			"MET CB", "MET CG", "MET SD", "MET CE", "ASN N", "ASN CA", "ASN C",
			"ASN O", "ASN CB", "ASN CG", "ASN OD1", "ASN ND2", "PRO N",
			"PRO CA", "PRO C", "PRO O", "PRO CB", "PRO CG", "PRO CD", "GLN N",
			"GLN CA", "GLN C", "GLN O", "GLN CB", "GLN CG", "GLN CD",
			"GLN OE1", "GLN NE2", "ARG N", "ARG CA", "ARG C", "ARG O",
			"ARG CB", "ARG CG", "ARG CD", "ARG NE", "ARG CZ", "ARG NH1",
			"ARG NH2", "SER N", "SER CA", "SER C", "SER O", "SER CB", "SER OG",
			"THR N", "THR CA", "THR C", "THR O", "THR CB", "THR OG1",
			"THR CG2", "VAL N", "VAL CA", "VAL C", "VAL O", "VAL CB",
			"VAL CG1", "VAL CG2", "TRP N", "TRP CA", "TRP C", "TRP O",
			"TRP CB", "TRP CG", "TRP CD1", "TRP CD2", "TRP NE1", "TRP CE2",
			"TRP CE3", "TRP CZ2", "TRP CZ3", "TRP CH2", "TYR N", "TYR CA",
			"TYR C", "TYR O", "TYR CB", "TYR CG", "TYR CD1", "TYR CD2",
			"TYR CE1", "TYR CE2", "TYR CZ", "TYR OH" };

	static double[][][] hphob_hphob_table = new double[atomType][atomType][maxCol];
	static double[][][] hphil_hphil_table = new double[atomType][atomType][maxCol];
	static double[][][] hphob_hphil_table = new double[atomType][atomType][maxCol];
	static double[][] freqTbl_hphob_hphob = new double[maxRow][maxCol];
	static double[][] freqTbl_hphil_hphil = new double[maxRow][maxCol];
	static double[][] freqTbl_hphob_hphil = new double[maxRow][maxCol];

	static double[][][] probTblHphob = new double[atomType][atomType][maxCol];
	static double[][][] probTblHphil = new double[atomType][atomType][maxCol];
	static double[][][] probTblHphob_Hphil = new double[atomType][atomType][maxCol];

	static int numAA = 20;
	static int numColASAEgyTable = 41;
	static double[][] ASA_Energy_Table = new double[numAA][numColASAEgyTable];
	static String[][] AA_Mapper = { { "A", "1" }, { "R", "2" }, { "N", "3" },
			{ "D", "4" }, { "C", "5" }, { "Q", "6" }, { "E", "7" },
			{ "G", "8" }, { "H", "9" }, { "I", "10" }, { "L", "11" },
			{ "K", "12" }, { "M", "13" }, { "F", "14" }, { "P", "15" },
			{ "S", "16" }, { "T", "17" }, { "W", "18" }, { "Y", "19" },
			{ "V", "20" } };

	static double predictedEnergy = 0.;
	
	
	// setup for uPhi and uPsi based energies
	static int maxDihedralCol = 20;		// bins range from -1 to 1 ... each bin of 0.1 radian ... total 20 bins
	static int maxCoord = 3;
	
	static double[][][] probTbl_dihedral_angle_phi = new double[atomType][atomType][maxDihedralCol];
	
	static double[][][] probTbl_dihedral_angle_psi = new double[atomType][atomType][maxDihedralCol];

	public static void main(String[] args) throws Exception {
		
		if(args.length == 0){
			
			System.out.println("Provide the valid pdbID. \nThe pdbID should match to the file name without extension to the file you placed in /Input/pdbInput/ directory within this software.\nE.g. valid pdbID is: 1ABBH");
			System.exit(0);
		}
		File currentDir = new File(new File(".").getAbsolutePath());
		String curr_dir_path = currentDir.getCanonicalPath();
		
		String log_file_writer = curr_dir_path+"/log_file.txt";
		Path log_path = Paths.get(log_file_writer);
		BufferedWriter log_writer = BufferReaderAndWriter.getWriter(log_path);
				
		
		String REGAd3pSoftwarePath = curr_dir_path+"/REGAd3p";
		String REGAd3pSoftwareScriptPath = REGAd3pSoftwarePath+"/Software/Scripts";
		String dsspPath = curr_dir_path;
		
		String pdbFileDir = curr_dir_path+"/Input/pdbInput/";
		
		String pdbFilePath = pdbFileDir+args[0]+".pdb";  // user needs to put pdb file in the pdbInput directory within the workspace and while running user is required to provide pdbID along with the java run command
	
		String pdbID = args[0];
		
		String firstRun = args[1];
		
		String[] fastaIdArray = pdbID.split("_");
			
		String fastaId = fastaIdArray[0];
		
		String rank = fastaIdArray[1];

		String consolidatedASA_Path = curr_dir_path+"/ConsolidatedASA/";

		int exitVal1 = runDSSPProcess(dsspPath, pdbFilePath, pdbID);

		if (exitVal1 == 0) {

			System.out.println("DSSP ran successfully");

		} else {

			System.err.println("DSSP failed to create DSSP file");

		}

		int exitVal2 = runDSSPOutputParser(curr_dir_path, pdbID, rank);
		if (exitVal2 == 0) {

			System.out.println("DSSP file parse successful");

		} else {

			System.err
					.println("DSSPOutputParser program failed to parse DSSP output");

		}
		
		if(firstRun.equals("T")){
			
			System.out.println("Running REGAd3p");
			int exitVal = runREGAd3pProcess(REGAd3pSoftwareScriptPath, rank);

			if (exitVal == 0) {

				System.out.println("REGAd3p run successful");

			} else {

				System.err.println("REGAd3p failed to produce predicted RSA");

			}
		
			System.out.println("Creating ASA file required for ASA Energy");
			createASAFileFor3DIGARS2(REGAd3pSoftwarePath, curr_dir_path, consolidatedASA_Path, pdbID, fastaId, rank, log_writer);

			System.out.println("ASA file required for ASA based energy is successfully created");
			
			
		}else if(firstRun.equals("F")){
			
			System.out.println("Creating ASA file required for ASA Energy");
			createASAFileFor3DIGARS2(REGAd3pSoftwarePath, curr_dir_path, consolidatedASA_Path, pdbID, fastaId, rank, log_writer);

			System.out.println("ASA file required for ASA based energy is successfully created");
			
		}
		

		readProbTableInMemory(curr_dir_path);

		System.out.println("3DIGARS library loaded successfully");

		readASAEnergyTable(curr_dir_path);
		System.out.println("ASA library loaded successfully");
		
		readuPhiProbTableInMemory(curr_dir_path);		
		System.out.println("uPhi Library Loaded Successfully");
		
		readuPsiProbTableInMemory(curr_dir_path);
		System.out.println("uPhi Library Loaded Successfully");
		
		
		System.out.println("Generating 3DIGARS3.0 Energy");
		generateThreeDIGARSVersionThreeEnergy(betas[0], betas[1], betas[2], w1, w2, w3,
				curr_dir_path, consolidatedASA_Path, pdbID, pdbFilePath, log_writer);

		// System.out.println("Done");
		log_writer.close();

	}

	public static void createASAFileFor3DIGARS2(String REGAd3pPath, String curr_dir_path, String consolidatedASA_Path,
			String pdbID, String fastaId, String rank, BufferedWriter log_writer) throws IOException {

		List<String> regad3pFileHolder = new ArrayList<>();
		List<String> dsspFileHolder = new ArrayList<>();
		String pred_rsa_path = REGAd3pPath+"/Software/Output/prediction/";
		String real_dsspOutput = curr_dir_path+"/FEATURES/";

		Path predicted_ASA_path = Paths.get(pred_rsa_path + fastaId + rank
				+ "/ASA", "/" + fastaId + ".ASAp");
		
		BufferedReader asa_pred_file_reader = BufferReaderAndWriter
				.getReaderNow(predicted_ASA_path);
				
		
		Path real_ASA_path = Paths.get(real_dsspOutput + pdbID + ".rsa");
		
		BufferedReader asa_real_file_reader = BufferReaderAndWriter
				.getReaderNow(real_ASA_path);
				
		

		String dsspLine = "";
		while ((dsspLine = asa_real_file_reader.readLine()) != null) {

			if (dsspLine.startsWith("#SR")) {

				String line = "";

				while ((line = asa_real_file_reader.readLine()) != null) {

					String[] tokens = line.split("\\s+");

					dsspFileHolder.add(tokens[0] + " " + tokens[1] + " "
							+ tokens[2]);

				}

			}

		}

		asa_real_file_reader.close();

		String regad3pLine = "";
		while ((regad3pLine = asa_pred_file_reader.readLine()) != null) {

			if (regad3pLine.startsWith("SR#")) {

				String line = "";

				while ((line = asa_pred_file_reader.readLine()) != null) {

					String[] tokens = line.split("\\s+");

					regad3pFileHolder.add(" " + tokens[2]);

				}

			}

		}

		asa_pred_file_reader.close();
		

		Path op_path = Paths.get(consolidatedASA_Path, pdbID
				+ ".rpASA");
		BufferedWriter op_writer = BufferReaderAndWriter.getWriter(op_path);
		op_writer.write("SN# AA ASAr ASAp");
		op_writer.newLine();
		
		
		
		
		
		if(regad3pFileHolder.size() != dsspFileHolder.size()){
			
			Path tempFilePath = Paths.get(curr_dir_path, "log_file_"+pdbID+"_"+fastaId+"_"+rank+".txt");
			BufferedWriter temp_writer = BufferReaderAndWriter.getWriter(tempFilePath);
			temp_writer.write("predicted ASA list Size " + regad3pFileHolder.size());
			temp_writer.newLine();
			temp_writer.write("real ASA list Size " + dsspFileHolder.size());
			temp_writer.newLine();
			String predicted_path = pred_rsa_path + fastaId + rank+ "/ASA"+"/" + fastaId + ".ASAp";
			temp_writer.write("predicted asa path "+predicted_path);
			temp_writer.newLine();
			String real_path = real_dsspOutput + pdbID + ".rsa";
			temp_writer.write("real asa path "+real_path);
			temp_writer.newLine();
			temp_writer.close();
			
		//	System.exit(1);
			
		}
		
		

		for (int i = 0; i < regad3pFileHolder.size(); i++) {

			op_writer.write(dsspFileHolder.get(i));
			op_writer.write(regad3pFileHolder.get(i));
			op_writer.newLine();

		}
		regad3pFileHolder.clear();
		dsspFileHolder.clear();
		op_writer.close();

	}

	static void generateThreeDIGARSVersionThreeEnergy(float b1, float b2,
			float b3, float w1, float w2, float w3, String curr_dir_path, String asaFilePath,
			String pdbID, String pdbFilePath, BufferedWriter log_writer) throws Exception {

		String asaFile = asaFilePath + pdbID + ".rpASA";
//		String inputDirPath = curr_dir_path + "/pdbInput";
//		File dir = new File(pdbFilePath);
//		File[] listOfFiles = dir.listFiles();
//		for (int i = 0; i < listOfFiles.length; i++) {
//
//			String filePath = listOfFiles[i].getAbsolutePath();
			findCAlphaAndLoadInMemory(pdbFilePath, log_writer, curr_dir_path);

			calcEucDistAndProtability(b1, b2, b3);

			Double ThreeDIGARSEgy = predictedEnergy / 100.;

			float EgyFromASA = calcASAEgyByFileName(asaFile);
			
			predictedEnergy = 0.;
			cAlphaCount = -1;
			findCAlphaAndLoadInMemoryForuPhiuPsi(pdbFilePath);
			calcDihedralAngleEnergyuPhi();
			double dihedral_energy_uphi = predictedEnergy;
			predictedEnergy = 0.;
			calcDihedralAngleEnergyuPsi();
			double dihedral_energy_upsi = predictedEnergy;

			Double combinedEnergy = ThreeDIGARSEgy + (w1 * EgyFromASA)+(w2 * dihedral_energy_uphi) + (w3 * dihedral_energy_upsi);
			System.out.println("# Individual_Energies "+ThreeDIGARSEgy+"\t"+(EgyFromASA)+"\t"+(dihedral_energy_uphi)+"\t"+(dihedral_energy_upsi));
			System.out.println("# "+pdbFilePath + ":-\t\t"
					+ combinedEnergy);
		

	}
	
public static void calcDihedralAngleEnergyuPhi() {
		
		for (int x = 0; x < cAlphaCount + 1; x++) {
			
			for (int y = x + 1; y < cAlphaCount + 1; y++) {

				if(y == cAlphaCount-1 || y == cAlphaCount-2) continue;		// ignore the last 2 atom in the structure
																
				
				if (resNum[x] == resNum[y])	continue;

				int col = 0;
				
				double eucDis = 0.;

				for (int i = 0; i < 3; i++) {

					double sqr = alphaCarbonArray[x][i]
							- alphaCarbonArray[y][i];

					eucDis += sqr * sqr;

				}

				eucDis = java.lang.Math.sqrt(eucDis);
				
				if(eucDis > 15) continue;
				
				double[] coord0 = new double[3];
				double[] coord1 = new double[3];
				double[] coord2 = new double[3];
				double[] coord3 = new double[3];
				double[] v1 = new double[3];
				double[] v2 = new double[3];
				double[] v3 = new double[3];
				double[] v4 = new double[3];
				double[] v5 = new double[3];
				float dotProd = 0f;
				double magV4 = 0.;
				double magV5 = 0.;
				double angle = 0.;
				
												
				/*
				 * The torsion angle way to get the atom for Dihedral angle calculation is
				 * get one atom from x and 3 atoms from y
				 */
				
				for(int i = 0; i < maxCoord; i++){		// obtain coordinates of 4 atoms to compute the uPhi
				
					coord0[i] = alphaCarbonArray[x][i];
					coord1[i] = alphaCarbonArray[y][i];
					coord2[i] = alphaCarbonArray[y+1][i];
					coord3[i] = alphaCarbonArray[y+2][i];
					
					v1[i] = coord1[i] - coord0[i];
					v2[i] = coord2[i] - coord1[i];
					v3[i] = coord3[i] - coord2[i];
				
				}
				
							
				v4 = crossProduct(v1, v2);
				v5 = crossProduct(v2, v3);
				
			
				dotProd = (float) dotProduct(v4, v5);
				
			
				magV4 = Math.sqrt(v4[0]*v4[0]+v4[1]*v4[1]+v4[2]*v4[2]);
				magV5 = Math.sqrt(v5[0]*v5[0]+v5[1]*v5[1]+v5[2]*v5[2]);
				
				if(magV4 == 0.0 || magV5 == 0.0){
					
					continue;
					
				}
				
				double angleComponent = dotProd/(magV4*magV5);

				DecimalFormat newFormat = new DecimalFormat("#.#####");
				double fiveDecimal =  Double.valueOf(newFormat.format(angleComponent));
				
				col = getBinIndex(fiveDecimal);
				
				if (col >= maxDihedralCol){
					System.out.println("Col "+col);
					System.exit(1);
				}				
				
				predictedEnergy += probTbl_dihedral_angle_phi[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];
				
			}

		}	
			
	}
	
	public static void calcDihedralAngleEnergyuPsi(){		

		
		for (int x = 0; x < cAlphaCount + 1; x++) {
			
			if(x == 0 || x == 1) continue;

			for (int y = x + 1; y < cAlphaCount + 1; y++) {

				if (resNum[x] == resNum[y])	continue;

				int col = 0;
				
				double eucDis = 0.;

				for (int i = 0; i < 3; i++) {

					double sqr = alphaCarbonArray[x][i]
							- alphaCarbonArray[y][i];

					eucDis += sqr * sqr;

				}

				eucDis = java.lang.Math.sqrt(eucDis);
				
				if(eucDis > 15) continue;
				
				double[] coord0 = new double[3];
				double[] coord1 = new double[3];
				double[] coord2 = new double[3];
				double[] coord3 = new double[3];
				double[] v1 = new double[3];
				double[] v2 = new double[3];
				double[] v3 = new double[3];
				double[] v4 = new double[3];
				double[] v5 = new double[3];
				float dotProd = 0f;
				double magV4 = 0.;
				double magV5 = 0.;
				double angle = 0.;
							
				
				/*
				 * The torsion angle way to get the atom for Dihedral angle calculation is
				 * get 3 atoms from x and 1 atoms from y
				 */
				
				for(int i = 0; i < maxCoord; i++){		// obtain coordinates of 4 atoms to compute the angle Psi
				
					coord0[i] = alphaCarbonArray[x-2][i];
					coord1[i] = alphaCarbonArray[x-1][i];
					coord2[i] = alphaCarbonArray[x][i];
					coord3[i] = alphaCarbonArray[y][i];
					
					v1[i] = coord1[i] - coord0[i];
					v2[i] = coord2[i] - coord1[i];
					v3[i] = coord3[i] - coord2[i];
				
				}
				
				v4 = crossProduct(v1, v2);
				v5 = crossProduct(v2, v3);
				
			
				dotProd = (float) dotProduct(v4, v5);

				
				magV4 = Math.sqrt(v4[0]*v4[0]+v4[1]*v4[1]+v4[2]*v4[2]);
				magV5 = Math.sqrt(v5[0]*v5[0]+v5[1]*v5[1]+v5[2]*v5[2]);
				
				if(magV4 == 0.0 || magV5 == 0.0){
					
					continue;
					
				}
				
				double angleComponent = dotProd/(magV4*magV5);
				
				DecimalFormat newFormat = new DecimalFormat("#.#####");
				double fiveDecimal =  Double.valueOf(newFormat.format(angleComponent));
				
				col = getBinIndex(fiveDecimal);
				
				if (col >= maxDihedralCol){
					System.out.println("Col "+col);
					System.exit(1);
				}				
				
				predictedEnergy += probTbl_dihedral_angle_psi[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];
			
			}

		}	
			
		
	}
	
	static int getBinIndex(double angle){
		int bin_index = -1;
		for(int i = 0; i < maxDihedralCol; i++){
			
			if(angle < -0.9 && angle >= -1.0){
				bin_index = 0;
			}else if(angle < -0.8 && angle >= -0.9){
				
				bin_index = 1;
				
			}else if(angle < -0.7 && angle >= -0.8){
				
				bin_index = 2;
				
			}else if(angle < -0.6 && angle >= -0.7){
				
				bin_index = 3;
				
			}else if(angle < -0.5 && angle >= -0.6){
				
				bin_index = 4;
				
			}else if(angle < -0.4 && angle >= -0.5){
				
				bin_index = 5;
				
			}else if(angle < -0.3 && angle >= -0.4){
				
				bin_index = 6;
				
			}else if(angle < -0.2 && angle >= -0.3){
				
				bin_index = 7;
				
			}else if(angle < -0.1 && angle >= -0.2){
				
				bin_index = 8;
				
			}else if(angle < 0.0 && angle >= -0.1){
				
				bin_index = 9;
				
			}else if(angle >= 0.0 && angle < 0.1){
				
				bin_index = 10;
				
			}else if(angle >= 0.1 && angle < 0.2){
				
				bin_index = 11;
				
			}else if(angle >= 0.2 && angle < 0.3){
				
				bin_index = 12;
				
			}else if(angle >= 0.3 && angle < 0.4){
				
				bin_index = 13;
				
			}else if(angle >= 0.4 && angle < 0.5){
				
				bin_index = 14;
				
			}else if(angle >= 0.5 && angle < 0.6){
				
				bin_index = 15;
				
			}else if(angle >= 0.6 && angle < 0.7){
				
				bin_index = 16;
				
			}else if(angle >= 0.7 && angle < 0.8){
				
				bin_index = 17;
				
			}else if(angle >= 0.8 && angle < 0.9){
				
				bin_index = 18;
				
			}else if(angle >= 0.9 && angle <= 1.0){
				
				bin_index = 19;
				
			}
			
		}
		
		return bin_index;
		
	}
	
	static double[] crossProduct(double[] v1, double[] v2){
		double[] vec = new double[3];
		
		double xcomp, ycomp, zcomp;
		
		xcomp = v1[1]*v2[2] - v1[2]*v2[1];
		ycomp = v1[2]*v2[0] - v1[0]*v2[2];
		zcomp = v1[0]*v2[1] - v1[1]*v2[0];
		
		vec[0] = xcomp;
		vec[1] = ycomp;
		vec[2] = zcomp;
		
		return vec;
		
		
	}
	
	static double dotProduct(double[] vA, double[] vB){
		
		double dotValue = 0.;
		
		for(int i=0; i < vA.length; i++){
			
			dotValue += vA[i]*vB[i]; 
		}
		
		return dotValue;
		
	}
	
	static void findCAlphaAndLoadInMemoryForuPhiuPsi(String fileName) throws IOException {

		File file = new File(fileName);	
	
//		System.out.println(fileName);
		BufferedReader fileReader = BufferReaderAndWriter.getReader(file);

		if(fileReader == null){
			
			System.out.println("error reading file:- "+fileName);
		}
		
		String atomName = null;
		String resName = null;

		String resCheck = "UNK";
		String distRes = null;
		int resID = -1;
		String ignore = null;

		String pdbFileLine = "";

		while ((pdbFileLine = fileReader.readLine()) != null) {

			if (!pdbFileLine.startsWith("ATOM")) {
				continue;
			}

			atomName = pdbFileLine.substring(13, 16).trim();
			if (atomName.startsWith("H"))
				continue;
			
//			if(atomName.startsWith("C")) continue; // igonre all carbon atoms

			resName = pdbFileLine.substring(17, 20).trim();

			ignore = resName + " " + atomName; // to check if residue name and
												// atom name pair is a valid
												// pair based on the resAtomPair
												// decleared above

			int found = 0;
			for (int i = 0; i < atomType; i++) {

				if (ignore.equals(resAtomPair[i])) {
					cAlphaCount++;
					found = 1;
					resAtomPairArrayID[cAlphaCount] = i;
					break;

				}

			}

			if (found == 0) {

				continue;
			}

			distRes = pdbFileLine.substring(17, 26).trim();

			if (!resCheck.equals(distRes)) {

				resCheck = distRes;

				resID++;

			}


			resNum[cAlphaCount] = resID;

		//	String xCor = pdbFileLine.substring(30, 38).trim();
		//	String yCor = pdbFileLine.substring(38, 46).trim();
		//	String zCor = pdbFileLine.substring(46, 54).trim();
			
		//	System.out.println(pdbFileLine);
		//	String xCor = pdbFileLine.substring(32, 39).trim();
		//	System.out.println(xCor);
		//	String yCor = pdbFileLine.substring(40, 48).trim();
		//	System.out.println(yCor);
		//	String zCor = pdbFileLine.substring(49, 57).trim();
		//	System.out.println(zCor);
		
		//	alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(xCor);
		//	alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(yCor);
		//	alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(zCor);
		
		//	String[] tokens = pdbFileLine.split("\\s+");
			
		//	alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(tokens[5]);
		//	alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(tokens[6]);
		//	alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(tokens[7]);

		
			String xCor = pdbFileLine.substring(30, 37).trim();
			String yCor = pdbFileLine.substring(38, 45).trim();
			String zCor = pdbFileLine.substring(46, 53).trim();
				
			alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(xCor);
			alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(yCor);
			alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(zCor);
					
		
		}

		fileReader.close();
		// printf("done reading file %s\n", fileName);
	}


	public static float calcASAEgyByFileName(String fileName)
			throws IOException {

		File file = new File(fileName);
		BufferedReader rosettaReader = BufferReaderAndWriter.getReader(file);
		String line = "";
		float ASA_Egy_ByFile = 0.f;
		while ((line = rosettaReader.readLine()) != null) {

			if (line.startsWith("SN#")) {
				continue;
			}

			String[] tokens = line.split("\\s+");

			String residue = tokens[1];

			int residueId = getAANumber(residue);
			if (residueId == 0) {

				continue;

			}
			float rsaR = Float.parseFloat(tokens[2]);

			float rsaP = Float.parseFloat(tokens[3]);

			float diff_rsaR_rsaP = java.lang.Math.abs(rsaR - rsaP);

			int bin = (int) java.lang.Math.floor(diff_rsaR_rsaP / 5);

			if (bin > 39) {

				bin = 39;
			}

			ASA_Egy_ByFile += ASA_Energy_Table[residueId - 1][bin + 1]; // array
																		// starts
																		// from
																		// 0,0
																		// and
																		// residueId
																		// is
																		// one
																		// ahead
																		// of
																		// array
																		// index,
																		// similarly
																		// bin
																		// starts
																		// from
																		// 1'st
																		// col,
																		// 0th
																		// col
																		// is AA
																		// index

		}

		return ASA_Egy_ByFile;

	}

	static void findCAlphaAndLoadInMemory(String fileName, BufferedWriter log_writer, String curr_dir_path) throws IOException, InterruptedException {

		try{
					
		File file = new File(fileName);
	//	log_writer.write(fileName);
	//	log_writer.write("\n");
		BufferedReader fileReader = BufferReaderAndWriter.getReader(file);
		
		if(fileReader == null){
			
			System.out.println(fileName);
		
		}

		String atomName = null;
		String resName = null;

		String resCheck = "UNK";
		String distRes = null;
		int resID = -1;
		String ignore = null;

		String pdbFileLine = "";

		while ((pdbFileLine = fileReader.readLine()) != null) {

			if (!pdbFileLine.startsWith("ATOM")) {
				continue;
			}

			atomName = pdbFileLine.substring(13, 16).trim();
			if (atomName.startsWith("H"))
				continue;

			resName = pdbFileLine.substring(17, 20).trim();

			ignore = resName + " " + atomName; // to check if residue name and
												// atom name pair is a valid
												// pair based on the resAtomPair
												// decleared above

			int found = 0;
			for (int i = 0; i < atomType; i++) {

				if (ignore.equals(resAtomPair[i])) {
					cAlphaCount++;
					found = 1;
					resAtomPairArrayID[cAlphaCount] = i;
					break;

				}

			}

			if (found == 0) {

				continue;
			}

			distRes = pdbFileLine.substring(17, 26).trim();

			if (!resCheck.equals(distRes)) {

				resCheck = distRes;

				resID++;

			}

			for (int i = 0; i < 20; i++) {

				if (Hphob_Hphil_Mapper[i][0].equals(resName)) {

					if (Hphob_Hphil_Mapper[i][1].equals("hphob")) {

						Hphob_Hphil_Array[cAlphaCount] = 1;

					} else if (Hphob_Hphil_Mapper[i][1].equals("hphil")) {

						Hphob_Hphil_Array[cAlphaCount] = 0;

					}

					break;

				}

			}

			resNum[cAlphaCount] = resID;
			
		//	log_writer.write(pdbFileLine);
		//	log_writer.write("\n");
			
		//	System.out.println(pdbFileLine);
			String xCor = pdbFileLine.substring(30, 37).trim();
			String yCor = pdbFileLine.substring(38, 45).trim();
			String zCor = pdbFileLine.substring(46, 53).trim();
		//	String[] tokens = pdbFileLine.split("\\s+");
			
		//	String xCor = pdbFileLine.substring(32, 39).trim();
		//	System.out.println(xCor);
		//	String yCor = pdbFileLine.substring(40, 48).trim();
		//	System.out.println(yCor);
		//	String zCor = pdbFileLine.substring(49, 57).trim();
		//	System.out.println(zCor);
		
		//	alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(tokens[5]);
		//	alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(tokens[6]);
		//	alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(tokens[7]);
		
			alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(xCor);
			alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(yCor);
			alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(zCor);

		}

		fileReader.close();
		
		}catch(NumberFormatException ex){
			
			copyThePdbWithError(fileName, curr_dir_path);
			
		}

	}
	
	public static void copyThePdbWithError(String fileName, String curr_dir_path) throws IOException, InterruptedException{
		
		String copy_cmd = "scp -r "+fileName+" "+curr_dir_path;
		Process p = Runtime.getRuntime().exec(copy_cmd);
		p.waitFor();		
		
	}

	public static void calcEucDistAndProtability(float b1, float b2, float b3) {

		predictedEnergy = 0.;
		double pred_hphob = 0.;
		double pred_hphil = 0.;
		double pred_hphob_hphil = 0.;

		for (int x = 0; x < cAlphaCount + 1; x++) {

			for (int y = x + 1; y < cAlphaCount + 1; y++) {

				if (resNum[x] == resNum[y])
					continue;

				int col = 0;

				double eucDis = 0.;

				for (int i = 0; i < 3; i++) {

					double sqr = alphaCarbonArray[x][i]
							- alphaCarbonArray[y][i];

					eucDis += sqr * sqr;

				}

				eucDis = java.lang.Math.sqrt(eucDis);
				col = (int) (eucDis * 2);

				if (col >= maxCol)
					continue;

				if (Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 1) {

					pred_hphob += probTblHphob[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				} else if (Hphob_Hphil_Array[x] == 0
						&& Hphob_Hphil_Array[y] == 0) {

					pred_hphil += probTblHphil[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				} else if ((Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 0)
						|| (Hphob_Hphil_Array[x] == 0 && Hphob_Hphil_Array[y] == 1)) {

					pred_hphob_hphil += probTblHphob_Hphil[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				}

			}

		}

		predictedEnergy = (b1 * pred_hphob) + (b2 * pred_hphil)
				+ (b3 * pred_hphob_hphil);

	}

	public static int runDSSPProcess(String dsspPath, String pdbFilePath,
			String pdbID) throws IOException, InterruptedException {
		System.out.println(pdbFilePath);
		ProcessBuilder pb1 = new ProcessBuilder("dssp", "-i", pdbFilePath);
		
		Map<String, String> env1 = pb1.environment();
		env1.put("VAR1", "myValue");
		env1.remove("OTHERVAR");
		env1.put("VAR2", env1.get("VAR1") + "suffix");
		pb1.directory(new File(dsspPath));
		File log1 = new File("./DSSP/" + pdbID + ".dssp");
		pb1.redirectErrorStream(true);
		pb1.redirectOutput(log1);
		Process p1 = pb1.start();
		assert pb1.redirectInput() == Redirect.PIPE;
		assert pb1.redirectOutput().file() == log1;
		assert p1.getInputStream().read() == -1;
		int exitVal1 = p1.waitFor();
		p1.destroy();
		return exitVal1;
		// System.out.println("Exited with error code "+exitVal1);

	}

	public static int runDSSPOutputParser(String curr_dir_path, String pdbID, String rank)
			throws IOException, InterruptedException {

		ProcessBuilder pb1 = new ProcessBuilder("./DSSPOutputParserMPI", pdbID);

		Map<String, String> env1 = pb1.environment();
		env1.put("VAR1", "myValue");
		env1.remove("OTHERVAR");
		env1.put("VAR2", env1.get("VAR1") + "suffix");
		pb1.directory(new File(curr_dir_path));
		File log1 = new File("logDSSPOutputParser_"+rank+".out");
		pb1.redirectErrorStream(true);
		pb1.redirectOutput(log1);
		Process p1 = pb1.start();
		assert pb1.redirectInput() == Redirect.PIPE;
		assert pb1.redirectOutput().file() == log1;
		assert p1.getInputStream().read() == -1;
		int exitVal1 = p1.waitFor();
		p1.destroy();
		return exitVal1;
		// System.out.println("Exited with error code "+exitVal1);

	}

	public static int runREGAd3pProcess(String REGAd3pSoftwareScriptPath, String rank)
			throws IOException, InterruptedException {

		ProcessBuilder pb = new ProcessBuilder("./run_REGAd3p", rank);
		Map<String, String> env = pb.environment();
		env.put("VAR1", "myValue");
		env.remove("OTHERVAR");
		env.put("VAR2", env.get("VAR1") + "suffix");
		pb.directory(new File(REGAd3pSoftwareScriptPath));
		File log = new File("logREGAd3p_"+rank+".out");
		pb.redirectErrorStream(true);
		pb.redirectOutput(log);
		Process p = pb.start();
		assert pb.redirectInput() == Redirect.PIPE;
		assert pb.redirectOutput().file() == log;
		assert p.getInputStream().read() == -1;
		int exitVal = p.waitFor();
		p.destroy();
		// System.out.println("Exited with error code "+exitVal);

		return exitVal;
	}

	public static void readProbTableInMemory(String curr_dir_path)
			throws Exception {

		String probTableHphob = curr_dir_path + "/probTableHphob.txt";
		String probTableHphil = curr_dir_path + "/probTableHphil.txt";
		String probTableHphob_Hphil = curr_dir_path
				+ "/probTableHphob_Hphil.txt";
		File file_hphob = new File(probTableHphob);
		BufferedReader br_hphob = BufferReaderAndWriter.getReader(file_hphob);
		File file_hphil = new File(probTableHphil);
		BufferedReader br_hphil = BufferReaderAndWriter.getReader(file_hphil);
		File file_hphob_hphil = new File(probTableHphob_Hphil);
		BufferedReader br_hphob_hphil = BufferReaderAndWriter
				.getReader(file_hphob_hphil);

		// load energy score in memory (probTableHphob.txt)
		String resAHphob;
		String resBHphob;

		String line_hphob = "";
		while ((line_hphob = br_hphob.readLine()) != null) {

			if (line_hphob.startsWith("#"))
				continue;

			String tokenHphob[] = line_hphob.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = tokenHphob[0];
			resBHphob = tokenHphob[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokenHphob.length; j++) {

				probTblHphob[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokenHphob[j]);
				probTblHphob[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokenHphob[j]);
			}

		}
		br_hphob.close();

		// load energy score in memory (probTableHphil.txt)

		resAHphob = null;
		resBHphob = null;

		String line_hphil = "";

		while ((line_hphil = br_hphil.readLine()) != null) {

			if (line_hphil.startsWith("#"))
				continue;

			String tokenHphil[] = line_hphil.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = tokenHphil[0];
			resBHphob = tokenHphil[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokenHphil.length; j++) {

				probTblHphil[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokenHphil[j]);
				probTblHphil[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokenHphil[j]);
			}

		}

		br_hphil.close();

		// load energy score in memory (probTableHphob_Hphil.txt)

		resAHphob = null;
		resBHphob = null;

		String line_hphob_hphil = "";

		while ((line_hphob_hphil = br_hphob_hphil.readLine()) != null) {

			if (line_hphob_hphil.startsWith("#"))
				continue;

			String token_Hphob_Hphil[] = line_hphob_hphil.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = token_Hphob_Hphil[0];
			resBHphob = token_Hphob_Hphil[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < token_Hphob_Hphil.length; j++) {

				probTblHphob_Hphil[resIDA][resIDB][j - 2] = Double
						.parseDouble(token_Hphob_Hphil[j]);
				probTblHphob_Hphil[resIDB][resIDA][j - 2] = Double
						.parseDouble(token_Hphob_Hphil[j]);
			}

		}

		br_hphob_hphil.close();

		// end of reading library file

	}

	public static int getAANumber(String st) {

		for (int i = 0; i < AA_Mapper.length; i++) {

			if (AA_Mapper[i][0].equalsIgnoreCase(st)) {

				return i + 1;

			}

		}

		return 0;

	}

	public static void readASAEnergyTable(String curr_dir_path)
			throws IOException {

		String ASAFilePath = curr_dir_path + "/EnergyTable_ln.csv";
		File ASA_Egy_File = new File(ASAFilePath);
		BufferedReader ASA_Egy_Reader = BufferReaderAndWriter
				.getReader(ASA_Egy_File);
		String line = "";
		int rowIndex = 0;
		while ((line = ASA_Egy_Reader.readLine()) != null) {

			String[] tokens = line.split("\\s+");

			ASA_Energy_Table[rowIndex][0] = getAANumber(tokens[0]);

			for (int j = 1; j < numColASAEgyTable; j++) {

				ASA_Energy_Table[rowIndex][j] = Double.parseDouble(tokens[j]);

			}

			rowIndex++;

		}

	}

	
public static void readuPhiProbTableInMemory(String curr_dir_path) throws Exception {
		
		String dihedralProbTablePath = curr_dir_path+"/probTableuPhi.txt";
		File file_dihedral = new File(dihedralProbTablePath);
		BufferedReader br_dihedral = BufferReaderAndWriter.getReader(file_dihedral);
		

		// load energy score in memory 
		
		String resA;
		String resB;
		
		String line_dihedral_energy = "";
		while ((line_dihedral_energy = br_dihedral.readLine()) != null) {

			if (line_dihedral_energy.startsWith("#"))
				continue;

			String tokens[] = line_dihedral_energy.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resA = tokens[0];
			resB = tokens[1];

			for (int i = 0; i < atomType; i++) {

				if (resA.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resB.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokens.length; j++) {

				probTbl_dihedral_angle_phi[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokens[j]);
				probTbl_dihedral_angle_phi[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokens[j]);
			}

		}
		br_dihedral.close();

		

	}
	
	
public static void readuPsiProbTableInMemory(String curr_dir_path) throws Exception {
		
		String dihedralProbTablePath = curr_dir_path+"/probTableuPsi.txt";
		File file_dihedral = new File(dihedralProbTablePath);
		BufferedReader br_dihedral = BufferReaderAndWriter.getReader(file_dihedral);
		

		// load energy score in memory 
		
		String resA;
		String resB;
		
		String line_dihedral_energy = "";
		while ((line_dihedral_energy = br_dihedral.readLine()) != null) {

			if (line_dihedral_energy.startsWith("#"))
				continue;

			String tokens[] = line_dihedral_energy.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resA = tokens[0];
			resB = tokens[1];

			for (int i = 0; i < atomType; i++) {

				if (resA.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resB.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokens.length; j++) {

				probTbl_dihedral_angle_psi[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokens[j]);
				probTbl_dihedral_angle_psi[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokens[j]);
			}

		}
		br_dihedral.close();

		

	}	
	

}
