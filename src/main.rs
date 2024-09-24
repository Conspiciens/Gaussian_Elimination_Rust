use std::fs::File; 
use std::io::{BufRead, BufReader}; 


fn fwd_elimination(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>) -> () {
    let _size = const_num.len(); 

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


fn back_subset(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, sol: &mut Vec<f64>) -> () {
    let _size = const_num.len(); 

    sol[_size - 1] = const_num[_size - 1] / coeff[_size - 1][_size - 1]; 

    for _i in (0.._size - 1).rev() {
        let mut sum: f64 = const_num[_i]; 
        for _j in _i + 1.._size { 
            sum = sum - coeff[_i][_j] * sol[_j]; 
        } 
        sol[_i] = sum / coeff[_i][_i]; 
    }         
}

fn naive_gaussian(mut coeff: Vec<Vec<f64>>, mut const_num: Vec<f64>) -> () {
    let size: usize = const_num.len();  
    let mut sol: Vec<f64> = vec![0.0; size]; 

    /* Borrowing the coeff and const_num */ 
    fwd_elimination(&mut coeff, &mut const_num);  
    back_subset(&mut coeff, &mut const_num, &mut sol); 
} 

/* Read from file the matrix and return the variable and const matrices */ 
fn load_matrix() -> (Vec<Vec<f64>>, Vec<f64>) {
    let file = File::open("sys1.lin").expect("File doesn't exist"); 
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
            .map(|num_str| num_str.trim().parse::<f64>().expect("Error parseing values"))
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

/* fn write_to_file(mut coeff: Vec<Vec<f64>>, mut const_num: Vec<f64>) -> None {
    let file = File::open(""); 

} */ 

fn spp_gausian(mut coeff: Vec<Vec<f64>>, mut const_num: Vec<f64>) -> () {
    let size: usize = coeff.len(); 
    let mut sol: Vec<f64>  = vec![0.0; size]; 
    let mut ind: Vec<usize> = vec![0; size]; 

    for _i in 0..coeff.len() {
        ind[_i] = _i; 
    } 

    sppd_fwd_elimination(&mut coeff, &mut const_num, &mut ind); 
    spp_back_subset(&mut coeff, &mut const_num, &mut sol, &mut ind); 
    
    println!("{:?}", sol); 

} 

fn sppd_fwd_elimination(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, ind: &mut Vec<usize>) -> () {
    let _size: usize = const_num.len(); 
    let mut scaling: Vec<f64> = vec![0.0; _size]; 
    
    for _i in 0.._size {
        let mut smax: f64 = 0.0; 
        for _j in 0.._size {
            if smax < coeff[_i][_j].abs() {
                smax = coeff[_i][_j].abs(); 
            } 
        }
        scaling[_i] = smax; 
    } 

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


fn spp_back_subset(coeff: &mut Vec<Vec<f64>>, const_num: &mut Vec<f64>, sol: &mut Vec<f64>, ind: &mut Vec<usize>) -> () {
    let _size: usize = const_num.len(); 

    sol[_size - 1] = const_num[ind[_size - 1]] / coeff[ind[_size - 1]][_size - 1]; 

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
    // let pattern = std::env::args().nth(1).expect("No flag was given"); 
    // let path = std::env::args().nth(2).expect("No File was given"); 


    /* Load the matrix from file */ 
    let matrix: (Vec<Vec<f64>>, Vec<f64>) = load_matrix(); 


    naive_gaussian(matrix.0, matrix.1); 
    spp_gausian(matrix.0, matrix.1); 
}
