//Angela Fan


var gene_name = "PRMT7";

//find the files the gene is in
d3.csv("/Documents/Bradner_work/hockey-sticks/lookup_table.csv", function(error, data){
	//MM1S_H3K27AC_DMSO_HOCKEY.csv 
	var input = gene_name;

	var data_length = data.length;

	var file_name_array = [];

	for (var i = 0; i < data.length; i++) {
		//console.log("hello")
		if (data[i].gene == input) {
			if (data[i].file) {
				file_name_array.push(data[i].file);
				console.log(data[i].file)
				break;
			}
		}
	}

	for (var i = 0; i < file_name_array.length; i++) {

		var current_file = file_name_array[i];

		d3.csv(current_file, function(error, current_data) {
			
			console.log(current_data)

			var current_file_name = current_file.splice("_HOCKEY")[0];

			current_file_name = "SE_plots_" + current_file_name;


		})

	}

});

// var pdf_name = "/Documents/Bradner_work/hockey-sticks/" + output_name[0] + "_plots/" + "SE_plots_" + output_name[0] + "_" + d.REGION_ID + ".pdf";

// $("#pdf_window")
// 	.append('<div class="pdf_image"> <object data=' + filename + 
// 		' type="application/pdf" width="800px" height="300"> alt : <a href='+ filename + '>test.pdf</a> </object> </div>' )                                                                               


