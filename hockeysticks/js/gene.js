//Angela Fan


var gene_name = "CKS1B";

//find the files the gene is in
d3.csv("/Documents/Bradner_work/hockey-sticks/lookup_table.csv", function(error, data){
	//MM1S_H3K27AC_DMSO_HOCKEY.csv 
	var input = gene_name;

	d3.select(".page_title").remove();

	d3.select("#page_title")
		.append("text")
		.attr("y", 0)
		.attr("text-anchor", "middle")
		.attr("x", 500)
		.text(gene_name)
		.attr("class", "page_title");

	var data_length = data.length;

	var file_name_array = [];

	for (var i = 0; i < data.length; i++) {
		//console.log("hello")
		if (data[i].gene == input) {
			if (data[i].file) {

				file_name_array.push(String(data[i].file).split(";"))

				//file_name_array.push(data[i].file);
				//console.log(data[i].file)
				break;
			}
		}
	}


	var ranking_array = [];

	for (var i = 0; i < file_name_array[0].length; i++) {

		var current_file = file_name_array[0][i];

		// console.log(i)

		// console.log(current_file)

		//console.log(current_file)

		plot(current_file, i)

		function plot(current_file, number) {
			d3.csv("/Documents/Bradner_work/hockey-sticks/" + current_file, function(error, current_data) {

				// console.log(current_file)

				// console.log(i)
				
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

							//console.log(current_file)


							ranking_array.push({"filename": current_file, "rank": +d.RANK, "prox_genes": d.PROXIMAL_GENES, 
								"super": +d.IS_SUPER, "signal": +d.SIGNAL, "top_gene": d.TOP_GENE, "function": d.PROXIMAL_FUNCTION,
								"start": +d.START, "stop": +d.STOP, "chromosome": d.CHROM});

							//console.log(ranking_array)

							//console.log(ranking_array)

							var current_file_name = current_file.split("_HOCKEY")[0];

							//console.log(d)

							current_file_name = "/Documents/Bradner_work/hockey-sticks/" + current_file_name + "_plots/" + "SE_plots_" + current_file_name + "_" + d.REGION_ID + ".pdf";

							$("#pdf_window")
								.append('<div class="gene_pdf"> <object data=' + current_file_name + 
								' type="application/pdf" width="840px" height="580"> alt : <a href='+ current_file_name + '>test.pdf</a> </object> </div>' )                                                                               


						}
					}

				});

				//console.log(ranking_array)

				// console.log(number)
				// console.log(file_name_array[0].length)

				if (number == 0) {

					tip = d3.tip().attr('class', 'd3-tip');

					var margin = {top: 50, right: 50, bottom: 50, left: 100},
					    width = 1100 - margin.left - margin.right,
					    height = 350 - margin.top - margin.bottom;

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

					var svg = d3.select("#bars").append("svg")
					    .attr("width", width + margin.left + margin.right)
					    .attr("height", height + margin.top + margin.bottom)
					  .append("g")
					    .attr("transform", "translate(" + (margin.left+100) + "," + margin.top+ ")");

					svg.call(tip)

					svg.append("text")
						.attr("x", 350)
						.attr("y", -20)
						.text("Comparison across all datasets")
						.attr("font-size", "18px")
						.attr("font-weight", "bold")

	
					function range(start, end) {
					    var foo = [];
					    for (var i = start; i <= end; i++) {
					        foo.push(i);
					    }
					    return foo;
					}


					x.domain(range(0,ranking_array.length-1));
					y.domain([0, d3.max(ranking_array, function(d) { return d.rank; })]);

					// console.log(x.domain())
					// console.log(x.range())

					svg.append("g")
				      	.attr("class", "x axis")
				      	.attr("transform", "translate(0," + height + ")")
				      	.call(xAxis)
				      	.append("text")
				      	.attr("x", 1000)
				      	.attr("y", 15)
				      	.style("text-anchor", "middle")
				      	.text("Files");

					svg.append("g")
				      	.attr("class", "y axis")
				      	.call(yAxis)
				    	.append("text")
				      	.attr("transform", "rotate(-90)")
				      	.attr("dy", "-3.71em")
				      	.style("text-anchor", "end")
				      	.text("Super-enhancer rank");

					svg.selectAll(".bar")
					    .data(ranking_array)
					    .enter().append("rect")
					    .attr("class", "bar")
					    .attr("x", function(d, j) {

					    	d.number = j
					    	//console.log(d)

					    	// console.log(j)
					    	// console.log(x(j))
					    	return x(j)
					    })
					    .attr("width", x.rangeBand())
					    .attr("y", function(d) { 

					    	return y(d.rank); 

					    })
					    .attr("height", function(d) { 

					    	return height - y(d.rank); 
					    })
					    .on("mouseover", function(d) {

					    	//console.log(d)

					    	tip.html("File: " + d.filename + "<br>Rank: " + d.rank + "<br>Signal: " + d.signal)

					    	return tip.show(d)

					    })
					    .on("mouseout", function(d) {
					    	return tip.hide(d)
					    });

					function type(d) {
					  d.frequency = +d.rank;
					  return d;
					}


					d3.select("input#by_rank").on("change", change_rank);
					d3.select("input#by_signal").on("change", change_signal);
					d3.select("input#by_super").on("change", change_super);

					function change_rank() {

						// Copy-on-write since tweens are evaluated after a delay.
						var x0 = x.domain(ranking_array.sort(this.checked
						    ? function(a, b) { return a.rank - b.rank; }
						    : function(a, b) { return d3.descending(a.number, b.number); })
						    .map(function(d) { 
						    	//console.log(d)
						    	return d.number; 
						    }))
						    .copy();

						var transition = svg.transition().duration(750),
						    delay = function(d, i) { return i * 50; };

						transition.selectAll(".bar")
						    .delay(delay)
						    .attr("x", function(d) { return x0(d.number); });

						transition.select(".x.axis")
						    .call(xAxis)
						  .selectAll("g")
						    .delay(delay);

					}

					function change_signal() {

						// Copy-on-write since tweens are evaluated after a delay.
						var x0 = x.domain(ranking_array.sort(this.checked
						    ? function(a, b) { return b.signal - a.signal; }
						    : function(a, b) { return d3.descending(a.number, b.number); })
						    .map(function(d) { 
						    	//console.log(d)
						    	return d.number; 
						    }))
						    .copy();

						var transition = svg.transition().duration(750),
						    delay = function(d, i) { return i * 50; };

						transition.selectAll(".bar")
						    .delay(delay)
						    .attr("x", function(d) { return x0(d.number); });

						transition.select(".x.axis")
						    .call(xAxis)
						  .selectAll("g")
						    .delay(delay);

					}

					function change_super() {

						// Copy-on-write since tweens are evaluated after a delay.
						var x0 = x.domain(ranking_array.sort(this.checked
						    ? function(a, b) { return a.super - b.super; }
						    : function(a, b) { return d3.descending(a.number, b.number); })
						    .map(function(d) { 
						    	//console.log(d)
						    	return d.number; 
						    }))
						    .copy();

						var transition = svg.transition().duration(750),
						    delay = function(d, i) { return i * 50; };

						transition.selectAll(".bar")
						    .delay(delay)
						    .attr("x", function(d) { return x0(d.number); });

						transition.select(".x.axis")
						    .call(xAxis)
						  .selectAll("g")
						    .delay(delay);

					}

				}

				
			})
		}


	}

});

