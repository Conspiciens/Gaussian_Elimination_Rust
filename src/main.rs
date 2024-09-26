use std::fs::File; 
use std::io::{ BufRead, BufReader, Write }; 
use std::env; 
/* 
    Within the file is code that will read a file that is specified by the CS3010 Assignment 1. Taking in a .lin file that contains the matrix size in row 1 for the columns and row     . The matrix numbers are spaced out equally and will be read by this program, the file must be in the same directory as the .rs file. If you want the binary file after running cargo run it's in /tar     /debug/gaussian. There are two flags that are avaliable to the user, which is --decimal to specify f64 or double percision or --spp to use the Sparse Partial Pivoting. Both can    can be used at the same time to activate both. The output will go to to .sol file and the solutions will be spaced out in a single line. The default is Gaussian elimination with single percision. : 
*/

/* 
    Write to file with the suffix being .sol 

    @params 
        file_name: Takes a filename input 
        sol: Solution Vec either in f32 or f64 (single or double per) 

    @return 
        None 
*/ 
fn write_to_file<T: ToString>(file_name: &String, sol: &Vec<T>) -> () {
    let mut file = File::create(file_name).expect("Error creating file"); 

    for num in sol {
        file.write_all(num.to_string().as_bytes()).expect("Error Writing to file"); 
        file.write_all(b" ").expect("Error Writing to file"); 
    } 
}  

/* 
    Foward Elimination using Naive Gaussian 
    
    @params
        coeff: Get all the coefficients 2d array 
        const_num: Get the const numbers 

    @return 
        None
*/ 
fn fwd_elimination(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>) -> () {
    let _size = const_num.len(); 

    /* Foward Elimination */ 
    for _k in 0.._size - 1 {    	
        for _i in _k + 1.._size {
            let multi: f64 = coeff[_i][_k] / coeff[_k][_k]; 
            for _j in _k.._size {
                coeff[_i][_j] -= multi * coeff[_k][_j]; 
            }  
            const_num[_i] = const_num[_i] - multi * const_num[_k]; 
        }
   }  
} 


/* 
   Back Substitution in the naive gaussian, get the final answers 
   @params 
        coeff: Coefficients stored within a 2d array 
        const_num: Constant numbers within a array 
        sol: Solutions derived from back substitution 

    @return 
        None
*/ 
fn back_subset(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, sol: &mut Vec<f64>) -> () {
    let _size = const_num.len(); 

    sol[_size - 1] = const_num[_size - 1] / coeff[_size - 1][_size - 1]; 

    /* Compute all the solutions */ 
    for _i in (0.._size - 1).rev() {
        let mut sum: f64 = const_num[_i]; 
        for _j in _i + 1.._size { 
            sum = sum - coeff[_i][_j] * sol[_j]; 
        } 
        sol[_i] = sum / coeff[_i][_i]; 
    }         
}

/* 
    Naive gaussian init the  
    @params
        coeff: Coefficients stored in 2d array
        const_num: Constant numbers stored in array   
    @return 
        Vec<f64>  
*/ 
fn naive_gaussian(mut coeff: Vec<Vec<f64>>, mut const_num: Vec<f64>) -> Vec<f64> {
    let size: usize = const_num.len();  
    let mut sol: Vec<f64> = vec![0.0; size]; 

    /* Borrowing the coeff and const_num */ 
    fwd_elimination(&mut coeff, &mut const_num);  
    back_subset(&mut coeff, &mut const_num, &mut sol); 

    return sol; 
} 

/* 
    Read from file the matrix and return the coefficients and const matrices 
    @params
        file_name: Name of the file to read from  
    @return 
        (Vec<Vec<f64>>, Vec<f64>)
*/ 
fn load_matrix(file_name: &String) -> (Vec<Vec<f64>>, Vec<f64>) {
    let file = File::open(file_name).expect("File doesn't exist"); 
    let reader = BufReader::new(file); 

    let mut total_lines: Option<usize> = None; 
    let mut matrix_vec: Vec<Vec<f64>> = Vec::new();  
    let mut const_vec: Vec<f64> = Vec::new(); 

    for (idx, line) in reader.lines().enumerate(){
    
        /* Split all the whitespace and convert to string */ 
        let vec_str: Vec<String> = line
            .expect("Line doesn't exist")
            .split_whitespace()
            .map(str::to_string)
            .collect(); 

        /* Convert all the strings to floats */ 
        let nums_array: Vec<f64> = vec_str 
            .into_iter()
            .map(|num_str| num_str.trim().parse::<f64>().expect("Error parsing values"))
            .collect(); 

        /* Init the matrix vector */ 
        if idx == 0 {
            total_lines = Some(nums_array[0] as usize); 
            matrix_vec = vec![vec![0.0; total_lines.expect("No usize")]; total_lines.expect("No usize")]; 
            continue; 
        } 

        /* Check lines is not None and is last line in file (const) to set const vec otherwise skip */ 
        match total_lines { 
            Some(line) if idx - 1 == line => {
                const_vec = vec![0.0; line]; 
                const_vec = nums_array; 
                break; 
            }, 
            _ => () 
        }; 
       
        /* Load the matrix with the correct nums */ 
        matrix_vec[idx - 1] = nums_array; 
    }    

    return (matrix_vec, const_vec); 
} 

/* 
    Scaled Partial Pivoting Gaussian Elimination Init 
    @params
        coeff: Coefficients in 2d array 
        const_num: constants within a array
    @return 
        sol: Vec<f64> 

*/ 
fn spp_gaussian(mut coeff: Vec<Vec<f64>>, mut const_num: Vec<f64>) -> Vec<f64> {
    let size: usize = coeff.len(); 
    let mut sol: Vec<f64>  = vec![0.0; size]; 
    let mut ind: Vec<usize> = vec![0; size]; 

    for _i in 0..coeff.len() {
        ind[_i] = _i; 
    } 

    spp_fwd_elimination(&mut coeff, &mut const_num, &mut ind); 
    spp_back_subset(&mut coeff, &mut const_num, &mut sol, &mut ind); 

    return sol; 
} 

/* 
    Scaled Partial Pivoting Gaussian Elimination Foward Elimination
    @params 
        coeff: Loads coefficients for 2d array 
        const_num: Loads the constant numbers in an array 
        ind: The indexes that are to be stored 

    @return 
        None 
*/ 
fn spp_fwd_elimination(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, ind: &mut Vec<usize>) -> () {
    let _size: usize = const_num.len(); 
    let mut scaling: Vec<f64> = vec![0.0; _size]; 
    
    /* Get the scaling numbers and find the greatest one */ 
    for _i in 0.._size {
        let mut smax: f64 = 0.0; 
        for _j in 0.._size {
            if smax < coeff[_i][_j].abs() {
                smax = coeff[_i][_j].abs(); 
            } 
        }
        scaling[_i] = smax; 
    } 

    /* Foward Elimination */ 
    for _k in 0.._size - 1 {
        let mut rmax: f64 = 0.0; 
        let mut maxInd: usize = _k; 

        for _i in _k.._size {
            let r: f64 = coeff[ind[_i]][_k] / scaling[ind[_i]]; 
            
            if r > rmax  {
                rmax = r; 
                maxInd = _i; 
            } 

        } 
    
        /* Swap */ 
        let tmp: usize = ind[maxInd]; 
        ind[maxInd] = ind[_k]; 
        ind[_k] = tmp;  

        for _i in _k + 1.._size {
            let mult: f64 = coeff[ind[_i]][_k] / coeff[ind[_k]][_k]; 
            
            for _j in _k + 1.._size { 
                coeff[ind[_i]][_j] = coeff[ind[_i]][_j] - mult * coeff[ind[_k]][_j];  
            } 

            const_num[ind[_i]] = const_num[ind[_i]] - mult * const_num[ind[_k]]; 
        } 
    } 

}  

/*
    Scaled Partial Pivoting Back Subtitution
    @params
        coeff: coefficients in a 2d array 
        const_num: constants that are within a arrray
        sol: solutions that are stored in an array
        ind: indexes that are used in switching

    @return 
        None
*/  
fn spp_back_subset(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, sol: &mut Vec<f64>, ind: &mut Vec<usize>) -> () {
    let _size: usize = const_num.len(); 

    sol[_size - 1] = const_num[ind[_size - 1]] / coeff[ind[_size - 1]][_size - 1]; 

    /* Compute all the solutions */ 
    for _i in (0..(_size - 1)).rev() {  
       let mut sum: f64 = const_num[ind[_i]]; 
       for _j in _i + 1.._size {
           sum = sum - coeff[ind[_i]][_j] * sol[_j];  
       }
       sol[_i] = sum / coeff[ind[_i]][_i]; 
    }  
}  

fn main() {
    /* Set the handling of spp flag and file input */ 
    let mut file_name: String= String::new();  
    let mut spp_flag: bool = false; 
    let mut double_flag: bool = false; 

    /* Collect arguments */ 
    let args: Vec<String> = env::args().collect(); 

    /* Check for arguments given by the user */ 
    for arg in args.into_iter() {
        match arg { 
            arg if arg.ends_with(".lin") => file_name = arg, 
            arg if arg.contains("--spp") => spp_flag = true, 
            arg if arg.contains("--double") => double_flag = true,  
            _ => () 
        } 
    } 

    /* Load matrix from file */ 
    let matrix: (Vec<Vec<f64>>, Vec<f64>) = load_matrix(&file_name); 

    /* Determine whether to run spp or naive */ 
    let mut sol: Vec<f64> = Vec::new(); 
    if spp_flag == true { 
        sol = spp_gaussian(matrix.0, matrix.1); 
    } else {
        sol = naive_gaussian(matrix.0, matrix.1); 
    } 

    /* Remove the .lin and add the .sol suffix */ 
    file_name = file_name.strip_suffix(".lin").expect("Error").to_string(); 
    file_name = file_name + ".sol"; 
    
    /* Check whether single or double percision is to be returned */ 
    if double_flag == true { 
        write_to_file(&file_name, &sol); 
    } else {
       let mut float_32_sol: Vec<f32> = Vec::new(); 
       for float_64 in &sol  {
            float_32_sol.push(*float_64 as f32); 
       } 
        write_to_file(&file_name, &float_32_sol); 
    } 
}
