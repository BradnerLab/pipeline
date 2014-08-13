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
				//console.log(data[i].file)
				break;
			}
		}
	}


	var ranking_array = [];

	for (var i = 0; i < file_name_array.length; i++) {

		var current_file = file_name_array[i];

		d3.csv("/Documents/Bradner_work/hockey-sticks/" + current_file, function(error, current_data) {
			
			current_data.forEach(function(d) {

				var genes = d.PROXIMAL_GENES.split(";");
				// if (i==1) {
				// 	console.log(genes)
				// }
				for (var j = 0; j < genes.length; j++) {
					if (i == 1) {
						//console.log(genes[j])
					}
					if (genes[j] == gene_name) {

						//console.log("hello")

						ranking_array.push({"filename": current_file, "rank": +d.RANK});

						//console.log(ranking_array)

						var current_file_name = current_file.split("_HOCKEY")[0];

						current_file_name = "/Documents/Bradner_work/hockey-sticks/" + current_file_name + "_plots/" + "SE_plots_" + current_file_name + "_" + d.REGION_ID + ".pdf";

						$("#pdf_window")
							.append('<div class="gene_pdf"> <object data=' + current_file_name + 
							' type="application/pdf" width="840px" height="580"> alt : <a href='+ current_file_name + '>test.pdf</a> </object> </div>' )                                                                               


					}
				}

			});

			//console.log(ranking_array)

			var margin = {top: 50, right: 50, bottom: 50, left: 100},
			    width = 1200 - margin.left - margin.right,
			    height = 300 - margin.top - margin.bottom;

			var x = d3.scale.ordinal()
			    .rangeRoundBands([0, width], .1);

			var y = d3.scale.linear()
			    .range([height, 0]);

			var xAxis = d3.svg.axis()
			    .scale(x)
			    .orient("bottom")
			    .tickFormat("");

			var yAxis = d3.svg.axis()
			    .scale(y)
			    .orient("left")
			    .outerTickSize([0]);

			var svg = d3.select("body").append("svg")
			    .attr("width", width + margin.left + margin.right)
			    .attr("height", height + margin.top + margin.bottom)
			  .append("g")
			    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

			x.domain([0, ranking_array.length]);
			y.domain([0, d3.max(ranking_array, function(d) { return d.rank; })]);

			svg.append("g")
		      	.attr("class", "x axis")
		      	.attr("transform", "translate(0," + height + ")")
		      	.call(xAxis);

			svg.append("g")
		      	.attr("class", "y axis")
		      	.call(yAxis)
		    	.append("text")
		      	.attr("transform", "rotate(-90)")
		      	.attr("y", 6)
		      	.attr("dy", ".71em")
		      	.style("text-anchor", "end")
		      	.text("Super-enhancer rank");

			svg.selectAll(".bar")
			    .data(ranking_array)
			    .enter().append("rect")
			    .attr("class", "bar")
			    .attr("x", function(d) { return x(d.filename); })
			    .attr("width", x.rangeBand())
			    .attr("y", function(d) { return y(d.rank); })
			    .attr("height", function(d) { return height - y(d.rank); });

			function type(d) {
			  d.frequency = +d.rank;
			  return d;
			}

		})

	}

});

