use std::fs::File; 
use std::io::{BufRead, BufReader}; 


fn fwd_elimination(coeff: &mut Vec<Vec<f32>>, const_num: &mut Vec<f32>){
    for _k in 0..coeff.len() - 1{    	
        let _i = _k + 1; 
        for _i in 0..coeff.len() {
            let multi = coeff[_i][_k] / coeff[_k][_k]; 
            let _j = _k; 
            for _j in 0..coeff.len() {
                coeff[_i][_j] = coeff[_i][_j] - multi * coeff[_k][_j]; 
            }  
            const_num[_i] = const_num[_i] - multi * const_num[_k]; 
        }
   }  
} 

fn back_subset(coeff: &mut Vec<Vec<f32>>, const_num: &mut Vec<f32>, sol: &mut Vec<f32>){
        
}

fn naive_gaussian(mut coeff: Vec<Vec<f32>>, mut const_num: Vec<f32>){
    let size: usize = coeff.capacity();  
    let mut sol: Vec<f32> = vec![0.0; size]; 
    
    /* Borrowing the coeff and const_num */ 
    fwd_elimination(&mut coeff, &mut const_num);  
    back_subset(&mut coeff, &mut const_num, &mut sol); 

} 

/* Read from file the matrix and return the variable and const matrices */ 
fn load_matrix() -> (Vec<Vec<f32>>, Vec<f32>) {
    let file = File::open("sys1.lin").expect("File doesn't exist"); 
    let reader = BufReader::new(file); 
    let mut total_lines: Option<usize> = None; 

    let mut matrix_vec: Vec<Vec<f32>> = Vec::new();  
    let mut const_vec: Vec<f32> = Vec::new(); 
    for (idx, line) in reader.lines().enumerate(){
    
        /* Split all the whitespace and convert to string */ 
        let vec_str: Vec<String> = line
            .expect("Line doesn't exist")
            .split_whitespace()
            .map(str::to_string)
            .collect(); 

        /* Convert all the strings to floats */ 
        let nums_array: Vec<f32> = vec_str 
            .into_iter()
            .map(|num_str| num_str.parse::<f32>().unwrap())
            .collect(); 


        /* Init the matrix vector */ 
        if idx == 0 {
            total_lines = Some(nums_array[0] as usize); 
            matrix_vec = vec![vec![0.0; total_lines.expect("No usize")]; total_lines.expect("No usize")]; 
            continue; 
        } 

        /* Check lines is not None and is last line to push const vec */ 
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

fn main() {
    let matrix: (Vec<Vec<f32>>, Vec<f32>) = load_matrix(); 

    naive_gaussian(matrix.0, matrix.1); 
}
