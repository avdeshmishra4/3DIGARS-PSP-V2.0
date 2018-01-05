//package FileOperations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;



public class BufferReaderAndWriter {
	
	// this method creates the buffer reader and return it back to the calling method
	public static BufferedReader getReaderNow(Path filePath) {

		Charset charset = Charset.forName("UTF-8");
		try {
			BufferedReader reader = Files.newBufferedReader(filePath, charset);
			return reader;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

	}
	
	public static BufferedReader getReader(File file) {
		FileInputStream fstream;
		try {
			FileInputStream fis = new FileInputStream(file);
			InputStreamReader isr = new InputStreamReader(fis, "UTF-8");
			BufferedReader br = new BufferedReader(isr);
			return br;

		} catch (Exception e) {
			return null;
		}
	}
	

	// this method creates the buffer writer and returns it back to the calling method
	
//	public static BufferedWriter getWriterNow(File file) {
//
//		FileOutputStream fstream;
//		try {
//			FileOutputStream fos = new FileOutputStream(file);
//			OutputStreamWriter osw = new OutputStreamWriter(fos, "UTF-8");
//			BufferedWriter bufferedWriter = new BufferedWriter(osw);
//			return bufferedWriter;
//
//		} catch (Exception e) {
//			return null;
//		}
//	}
	
	public static BufferedWriter getWriter(Path filePath) {
		Charset charset = Charset.forName("UTF-8");
		
		try {
			BufferedWriter writer = Files.newBufferedWriter(filePath, charset);
			return writer;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		
	}
	

}
