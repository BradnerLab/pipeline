//Angela Fan

//global variable for data
var data;

var file;

//margins and bounding boxes for the graph
var bb_linegraph;

bb_linegraph = {
	h: 300,
	w: 720,
	margin: {
		top: 50,
		right: 50,
		bottom: 20,
		left: 50
	}
};

$(document).ready(function() {
    $('[title]').removeAttr('title');
});

$(".dropdown-menu li a").click(function () {

	d3.select(".bubble").remove();
	d3.select(".crc_svg").remove();		
	d3.select("#slidertext").text("");
	d3.select("#slider").html("");
	d3.select(".linegraph_svg").remove();
	d3.select("#table_div")
		.style("visibility", "visible");

	d3.select("#add_all_visible")
		.style("visibility", "visible");

	d3.select(".linegraph_clear_button").remove();

	file = $(this).text()
	update_linegraph(file);

	output_name = file.split("_HOCKEY");

	d3.select(".page_title").remove();

	d3.select("#page_title")
		.append("text")
		.attr("y", 0)
		.attr("text-anchor", "middle")
		.attr("x", 500)
		.text(file.split("_HOCKEY")[0])
		.attr("class", "page_title");
});

//tip call
var graph_tip = d3.tip()
	.attr("class", "d3-tip")
	.offset([0,0]);

//linegraph globals
var linegraph_x, linegraph_y, linegraph_color, linegraph_xAxis, linegraph_yAxis;
var linegraph_xdomain, linegraph_ydomain;

function draw_linegraph() {

	svg_linegraph = d3.select("#linegraph").append("svg")
	.attr("class", "linegraph_svg")
	.attr({
		width: bb_linegraph.w + bb_linegraph.margin.left + bb_linegraph.margin.right + 400,
		height: bb_linegraph.h + bb_linegraph.margin.bottom + bb_linegraph.margin.top
	});

	svg_linegraph.call(graph_tip);

	var linegraph = svg_linegraph.append("g")
		.attr("class", "linegraph")
		.attr("transform", "translate(" + (bb_linegraph.margin.left-40) + "," + bb_linegraph.margin.top + ")");

	linegraph_x = d3.scale.linear()
		.range([bb_linegraph.w, 0]);

	linegraph_y = d3.scale.linear()
		.range([bb_linegraph.h, 0]);

	linegraph_color = d3.scale.ordinal()
		.domain([0,1])
		.range(["#868686", "#b30000"]);

	linegraph_xAxis = d3.svg.axis()
		.scale(linegraph_x)
		.orient("bottom")
		//get rid of weird default end tick
		.outerTickSize([0]);

	linegraph_yAxis = d3.svg.axis()
		.scale(linegraph_y)
		.orient("right")
		//get rid of weird default end tick
		.outerTickSize([0]);

	//add the clip path for the brush
	linegraph.append("defs")
		.append("clipPath")
		.attr("transform", "translate(-4.5,-4.5)")
		.attr("id", "linegraph_clip")
		.append("rect")
		.attr("width", bb_linegraph.w+4.5)
		.attr("height", bb_linegraph.h + 4.5);

	linegraph.append("g")
		.attr("class", "linegraph_brush");

	linegraph.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + bb_linegraph.h + ")")
		.call(linegraph_xAxis)
		.append("text")
		.attr("class", "label")
		.attr("y", 38)
		.attr("x", bb_linegraph.w/2 + 50)
		.attr("dy", "-.71em")
		.style("text-anchor", "end")
		.text("Enhancers Ranked by Signal");

	linegraph.append("g")
		.attr("class", "y axis")
		.attr("transform", "translate(" + bb_linegraph.w + ",0)")
		.call(linegraph_yAxis)
		.append("text")
		.attr("class", "label")
		.attr("transform", "rotate(-90)")
		.attr("y", 6)
		.attr("dy", "6.5em")
		.attr("x", -115)
		.style("text-anchor", "end")
		.text("Enhancer Signal");

	//titles 
	linegraph.append("text")
		.attr("y", -38)
		.attr("font-size", "16px")
		.attr("font-weight", "bold")
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.text("Hockey Stick");

	linegraph.append("text")
		.attr("y", -21)
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.style("font-size", "12px")
		.text("Select a region to zoom in. Click 'Clear Zoom' to zoom out. Click the buttons to add points to the table.");

	linegraph.append("text")
		.attr("y", -8)
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.style("font-size", "12px")
		.text("Mouseover a circle more information, click to bring up a PDF, double-click to add to the table.");

	d3.csv("/Documents/Bradner_work/hockey-sticks/lookup_table.csv", function(error, data){

		$("#search_button").click(function() {
			//console.log("hello")
			var input = $("#search").val()

			var data_length = data.length;

			for (var i = 0; i < data.length; i++) {
				if (data[i].gene == input) {
					if (data[i].file) {
						alert(data[i].file)
						break;
					}
					else {
						alert("Sorry " + input + " is in the dataset but is not super.")
						break;
					}
				}
			}

		})



	})

};

d3.selection.prototype.moveToFront = function() {
  return this.each(function(){
    this.parentNode.appendChild(this);
  });
};


var linegraph_brush, linegraph_xdomain, linegraph_ydomain;

function update_linegraph(file) {

	draw_linegraph();

	$("#exit").on("click", function() {
		d3.select("#linegraph").html("");
		clearTable("tbody");		
		d3.select("#table_div")
			.style("visibility", "hidden");
		d3.select("#crcgraph").html("");
		d3.select("#functional_bubble").html("");
		d3.select("#page_title").html("");
		d3.select("#add_all_visible")
			.style("visibility", "hidden");
	})

	$("#minimize").click(function(){
	    if($(this).html() == "Minimize"){
	        $(this).html("Maximize");
	    }
	    else{
	        $(this).html("Minimize");
	    }
	    $("#dataset_total").slideToggle();
	});

	$(".tablesorter").trigger("update");

	//reset the clear button 
	d3.select(".clear_button_linegraph").remove();

	line = d3.svg.line()
		.x(function(d) {return linegraph_x(d.rank)})
		.y(function(d) {return linegraph_y(d.val)});

	//data
	d3.csv("/Documents/Bradner_work/hockey-sticks/" + file, function(error, data) {
		data.forEach(function(d) {
			d.super = +d.IS_SUPER;
			d.val = +d.SIGNAL;
			d.rank = +d.RANK;
		});

		var data_length = data.length;

		make_graphs(data);

	    d3.select('#slider').call(d3.slider().axis(false).min(1).max(data_length).step(1)
	    	.on("slide", function(evt, value) {

	    		if (value != data_length) {

		    		d3.csv("/Documents/Bradner_work/hockey-sticks/" + file, function(error, data) {
						data.forEach(function(d) {
							d.super = +d.IS_SUPER;
							d.val = +d.SIGNAL;
							d.rank = +d.RANK;
						});

		      			d3.select('#slidertext').text(data_length - value);

		      			data.splice(data_length - value, value);

		      			make_graphs(data);

					})
	    		}
	    	})
	    );

		function make_graphs(data) {

      		d3.selectAll(".circle").remove();
      		d3.selectAll(".legend").remove();
      		d3.selectAll(".bubble").remove();

			var possible_categories = ["BloodSurfaceAntigen", "Bromodomain", "Cancer_Mutated", "Cardio_disease", "CD_Marker", "Chromatin_Modifying", "GPCR", "Kinase", "Methyltransferase", "Peptidases", "PHD_Containing", "Ribosomal", "Secreted_Hormone", "Transporter", "Tudo_Containing", "TxnFactor"];

			var functional_dict = [];

			var number_categories = possible_categories.length;


			for (var i = 0; i < number_categories; i++) {
				functional_dict.push({
					name: possible_categories[i],
					size: 0
				});
			}

			//console.log(data)
			//console.log(data)

			for (var i = 0; i < data.length; i++) {


				var current_categories = data[i].PROXIMAL_FUNCTION
				//console.log(data[i])
				var split_categories = current_categories.split(";");

				// if (!current_categories) {
				// 	console.log(data[i])
				// }

				var num_split = split_categories.length;

				for (var j = 0; j < num_split; j++) {

					for (var k = 0; k < number_categories; k++) {
						if (split_categories[j] == functional_dict[k].name) {
							functional_dict[k].size += 1;
						}

					}
				}
			}

			var flare_dict = {
				name: "flare",
				children: functional_dict
			};

			//console.log(flare_dict)

			var flare_diameter = 300,
			    flare_format = d3.format(",d");

			var flare_color = d3.scale.ordinal()
			    .domain(d3.extent(functional_dict, function(d) {
			    	return d.size;
			    }))
			    .range(["#3288bd","#66c2a5","#abdda4","#e6f598","#fee08b","#fdae61","#f46d43","#d53e4f"]);

			var bubble = d3.layout.pack()
			    .sort(null)
			    .size([flare_diameter, flare_diameter])
			    .padding(1.5);

			var functional_svg = d3.select("#functional_bubble").append("svg")
			    .attr("width", flare_diameter+80)
			    .attr("height", flare_diameter+80)
			    .attr("class", "bubble");

			var func_tip = d3.tip()
				.attr("class", "d3-tip")
				.offset([0,0]);

			functional_svg.call(func_tip);

			//console.log(bubble.nodes(classes(flare_dict)))

			var functional_node = functional_svg.selectAll(".flare_node")
			      .data(bubble.nodes(classes(flare_dict))
			      .filter(function(d) { return !d.children; }))
			    .enter().append("g")
			      .attr("class", "flare_node")
			      .attr("transform", function(d) { 
			      	//console.log(d)
			      	return "translate(" + d.x + "," + (d.y+80) + ")"; });

			functional_node.append("title")
			      .text(function(d) { return d.className + ": " + flare_format(d.value); });

			d3.selectAll("title").remove()

			functional_node.append("circle")
			      .attr("r", function(d) { return d.r; })
			      .style("fill", function(d) { return flare_color(d.className); })
			      .on("click", function(d) {
			      	for (var i = 0; i < data.length; i++) {

							current_categories = data[i].PROXIMAL_FUNCTION

							var split_categories = current_categories.split(";");

							var num_split = split_categories.length;

							for (var j = 0; j < num_split; j++) {

								for (var k = 0; k < number_categories; k++) {
									if (split_categories[j] == d.className) {
										//console.log("here")
										addRow("tbody", data[i])
										break;

									}

								}
							}
						}
			      })
			      .on("mouseover", function(d) {

			      	func_tip.html(d.className + ", " + d.value)
			      	func_tip.show(d);

						for (var i = 0; i < data.length; i++) {

							current_categories = data[i].PROXIMAL_FUNCTION

							var split_categories = current_categories.split(";");

							var num_split = split_categories.length;

							for (var j = 0; j < num_split; j++) {

								for (var k = 0; k < number_categories; k++) {
									if (split_categories[j] == d.className) {
										//console.log("here")

										d3.select("[rank='" + data[i].RANK + "']")
											.moveToFront()
											.attr("r", 6.5)
											.attr("stroke", "black")
											.attr("fill", "#354299")
											.attr("stroke-width", "2px")

									}

								}
							}
						}

			      })
			      .on("mouseout", function(d) {
			      	func_tip.hide(d);

			      	for (var i = 0; i < data.length; i++) {

							current_categories = data[i].PROXIMAL_FUNCTION

							var split_categories = current_categories.split(";");

							var num_split = split_categories.length;

							for (var j = 0; j < num_split; j++) {

								for (var k = 0; k < number_categories; k++) {
									if (split_categories[j] == d.className) {
										//console.log("here")

										d3.select("[rank='" + data[i].RANK + "']")
											.attr("stroke", null)
											.attr("fill", linegraph_color(data[i].super))
											.attr("r", 4.5);

									}

								}
							}
						}
			      });

			functional_node.append("text")
			      .attr("dy", ".3em")
			      .style("text-anchor", "middle")
			      .style("font-weight", "bold")
			      .style("font-size", "12px")
			      .text(function(d) { 
			      	//console.log(d.className.length, d.r/3)
			      	if (d.className.length < d.r/3) {
			      		return d.className;
			      	}
			      	//return d.className.substring(0, d.r / 3); 
			      });


			functional_svg.append("text")
				.attr("y", 40)
				.attr("text-anchor", "middle")
				.attr("x", 150)
				.style("font-size", "12px")
				.text("Mouseover for more information.");

			functional_svg.append("text")
				.attr("y", 55)
				.attr("text-anchor", "middle")
				.attr("x", 150)
				.style("font-size", "12px")
				.text("Click to add all genes with that annotation to the table");
			 //Click to add all genes with that annotation to the table.");


			function classes(root) {

			  var classes = [];

			  function recurse(name, node) {
			    if (node.children) node.children.forEach(function(child) { recurse(node.name, child); });
			    else classes.push({packageName: name, className: node.name, value: node.size});
			  }

			  recurse(null, root);
			  return {children: classes};

			}

			functional_svg.append("text")
				.attr("x", 150)
				.attr("y", 20)
				.attr("font-size", "16px")
				.attr("font-weight", "bold")
				.attr("text-anchor", "middle")
				.text("Functional Annotations");

			//console.log(data.length)

			linegraph_x.domain([0, data.length]);

			linegraph_y.domain(d3.extent(data, function(d) {
				return d.val;
			}));


			var linegraph = d3.select(".linegraph");

			var datapoints = linegraph.selectAll(".dot")
				.data(data, function(d) {return d.rank});

			datapoints.exit().remove();


			//add rows to table
			function addRow(tableID, d) {

				var table = document.getElementById(tableID);

				//check if the dot has been added before. 
				var table_row = $("#table tbody tr");

				var rank_array = [];

				table_row.each(function(index, element) {
					var rank_row = $(this).find("td").eq(2).text();
					rank_array.push(rank_row);
				});

				//if not, add it
				if ($.inArray(d.rank, rank_array) == -1) {
					// console.log(table)

					var rowCount = table.rows.length;
					var row = table.insertRow(rowCount);

					var cell0 = row.insertCell(0);
					var element1 = document.createElement("input");
					element1.type = "checkbox";
					element1.name = "check";
					cell0.appendChild(element1);

					var celladd = row.insertCell(1)
					celladd.innerHTML = output_name[0];

					var cell1 = row.insertCell(2);
					var split = d.REGION_ID.split("_");
					if (split[3] == "lociStitched") {
						cell1.innerHTML = split[0];
					}
					else {
						cell1.innerHTML = split[3];
					}

					// console.log(d.peak)

					var cell2 = row.insertCell(3);
					cell2.innerHTML = d.RANK;

					var cell3 = row.insertCell(4);
					if (d.super == 1) {
						cell3.innerHTML = "Yes";
						cell3.style.backgroundColor = "#b30000";
					}
					else {
						cell3.innerHTML = "No";
					}

					var cell4 = row.insertCell(5);
					cell4.innerHTML = d.CHROM;

					var cell5 = row.insertCell(6);
					cell5.innerHTML = d.START;

					var cell6 = row.insertCell(7);
					cell6.innerHTML = d.STOP;

					var cell7 = row.insertCell(8);
					cell7.innerHTML = d.PROXIMAL_GENES;

					var cell8 = row.insertCell(9);
					cell8.innerHTML = d.PROXIMAL_FUNCTION;

					//highlight rows on the table on mouseover
					var rows = d3.selectAll("tbody tr")
						.on("mouseover", function() {
							d3.select(this).style("background-color", "#354299");

							var row_rank = $(this).find("td").eq(3).text();

							var here = d3.select("[rank='" + row_rank + "']")
								.moveToFront()
								.attr("r", 6.5)
								.attr("stroke", "black")
								.attr("fill", "#354299")
								.attr("stroke-width", "2px");

						})
						.on("mouseout", function() {
							d3.select(this).style("background-color", null);

							var row_rank = $(this).find("td").eq(3).text();

							var row_super = $(this).find("td").eq(4).text();

							if (row_super == "No") {
								var superness = 0;
							}
							else {
								var superness = 1;
							}

							var here = d3.select("[rank='" + row_rank + "']")
								.attr("stroke", null)
								.attr("fill", linegraph_color(superness))
								.attr("r", 4.5);
						});

					//sorting the table
					$(function() {
						$('#table').tablesorter({
							headers: {
								0: {
									sorter: false
								},
								1: {
									sorter: false
								},
								2: {
									sorter: false
								},
								5: {
									sorter: false
								},
								6: {
									sorter: false
								},
								7: {
									sorter: false
								},
								8: {
									sorter: false
								},
								9: {
									sorter: false
								}
							}
						})
					})
				};  
			}

			linegraph.select(".x.axis")
				.transition()
				.duration(1000)
				.call(linegraph_xAxis);

			linegraph.select(".y.axis")
				.transition()
				.duration(1000)
				.call(linegraph_yAxis);


			var delay = 700, clicks = 0, timer = null;

			datapoints
				.enter()
				.append("circle")
				.attr("class", "dot")
				.attr("rank", function(d) {
					return d.rank;
				})
				.attr("super", function(d) {
					return d.super;
				})
				.attr("r", 4.5)
				.attr("opacity", "0.9")
				.attr("clip-path", "url(#linegraph_clip)")
				.on("mouseover", function(d) {

					//console.log(d)

					d3.select(this)
					.moveToFront()
					.attr("r", 6.5)
					.attr("stroke", "black")
					.attr("stroke-width", "2px");

					if (!d.PROXIMAL_GENES) {
						var nearby = "None";
					}
					else {
						var nearby = d.PROXIMAL_GENES;
					}
					//console.log(d)

					if (d[""]) {
						var gene = d[""];
					}
					else if(!d.top) {
						var gene = "None";
					}
					else {
						var gene = d.top;
					}

					// console.log(d[""])

					var split = d.REGION_ID.split("_");
					//console.log(d)
					if (split[3] == "lociStitched") {
						var peak = split[0];
					}
					else {
						var peak = split[3];
					}

					graph_tip.html("<strong>Peak Number: </strong>" + peak + "<br><strong>Top Gene: </strong>" + gene + "<br><strong>Chromosome: </strong>" + d.CHROM + "<br><strong>Nearby genes: </strong>" + nearby);
					graph_tip.show(d);
				})
				.on("mouseout", function(d) {
					d3.select(this)
					.attr("r", 4.5)
					.attr("stroke", null)

					graph_tip.hide(d);
				})
				.on("click", function(d) {
					clicks++;

					//console.log(d)

					if (clicks === 1 & d.super == 1) {

						//console.log(file)
	//output_name[0] + "plots/" +
						var pdf_name = "/Documents/Bradner_work/hockey-sticks/" + output_name[0] + "_plots/" + "SE_plots_" + output_name[0] + "_" + d.REGION_ID + ".pdf";

						//console.log(output_name)

						//console.log(pdf_name)
						//U87_H3K27AC_plots/SE_plots_U87_H3K27AC_1_MACS_peak_6695_lociStitched.pdf
						//U87_H3K27AC_plots/SE_plots_U87_H3K27AC6_MACS_peak_4435_lociStitched.pdf

						timer = setTimeout(function() {
							clicks = 0;
							//console.log("happen")
							screenshotPreview(pdf_name);

							// d3.selectAll(".pdf_image", function(d) {
							// 	console.log("here")
							// 	this.attr("transform", "translate(" + d.rank + "0)");
							// })

						}, delay);
					}

					else if (clicks === 1) {
						timer = setTimeout(function() {
							clicks = 0;
						}, delay);
					}

					else {

						clearTimeout(timer);
						clicks = 0
						addRow("tbody", d);
					}
				});

				d3.select("#linegraph")
					.on("click", function() {
						d3.select(".pdf_image").remove();
					})

			// console.log(data)

			datapoints
				.attr("cx", linegraph_x(1))
				.attr("cy", linegraph_y(1))
				.transition()
				.duration(1000)
				.attr("cx", function(d) {
					return linegraph_x(d.rank);
				})
				.attr("cy", function(d) {
					return linegraph_y(d.val);
				})
				.attr("fill", function(d) {
					return linegraph_color(d.super);
				});


			//set up the color legend
			var legend = linegraph.append("g")
				.attr("class", "legend")
				.attr("width", 200)
				.attr("height", 200)
				.selectAll("g")
				.data(linegraph_color.domain().slice().reverse())
				.enter()
				.append("g")
				.attr("transform", function(d,i) {
					return "translate(-10," + -1.5*(i*20) + ")";
				});

			legend.append("rect")
				.attr("width", 18)
				.attr("height", 18)
				.style("fill", linegraph_color);

			legend.append("text")
				.attr("x", 25)
				.attr("y", 0)
				.attr("dy", "1.15em")
				.text(function(d) {
					if (d == 1) {
						return "Super-enhancer"
					}
					else {
						return "Typical Enhancer"
					}
				});

			linegraph.append("text")
				.attr("x", 25)
				.attr("y", -40)
				.text("Legend")
				.attr("class", "legend")
				.style("font-weight", "bold")
				.style("font-size", "14px");

			//set up the brush
			linegraph_brush = d3.svg.brush()
				.x(linegraph_x)
				.y(linegraph_y)
				.on("brushend", linegraph_brushend);

			var max_x = d3.max(data, function(d) {
				return d.rank;
			});
			var max_y = d3.max(data, function(d) {
				return d.val;
			});

			//console.log(max_y)

			linegraph_xdomain = [0, max_x];
			linegraph_ydomain = [0, max_y];

			linegraph_x.domain(linegraph_xdomain);
			linegraph_y.domain(linegraph_ydomain);

			d3.select(".linegraph_brush").call(linegraph_brush);

			//updating the brush based on selected data

			var linegraph_clear_button;

			function linegraph_brushend() {

				var linegraph_x_domain = [linegraph_brush.extent()[0][0], linegraph_brush.extent()[1][0]]
				var linegraph_y_domain = [linegraph_brush.extent()[0][1], linegraph_brush.extent()[1][1]]

				// console.log(linegraph_x_domain)
				// console.log(linegraph_y_domain)

				// equal domain ends means click on graph
				if (+linegraph_x_domain[0] == +linegraph_x_domain[1] && +linegraph_y_domain[0] == +linegraph_y_domain[1]) {
					return;
				}

				get_button = d3.select(".linegraph_clear_button");

				if (get_button.empty() === true)
				{
					linegraph_clear_button = linegraph.append("g")
						.attr("transform", "translate(" + (bb_linegraph.w - 100) + "," + (bb_linegraph.h - 320) + ")")
						.attr("class", "linegraph_clear_button");

					linegraph_clear_button.append("rect")
						.attr("width", 102)
						.attr("height", 20)
						.attr("y", -20)
						.attr("x", -4)
						.attr("rx", "10px")
						.attr("ry", "10px")
						.style("fill", "#9f9f9f");

					linegraph_clear_button
						.append('text')
						.attr("y", -3)
						.attr("x", 10)
						.text("Clear Zoom")
						.style("fill", "black");
				}

				linegraph_x.domain(linegraph_x_domain);
				linegraph_y.domain(linegraph_y_domain);

				//add_to_table();

				$("#add_all_visible").click(function () {

					return add_to_table();
				});

				function add_to_table() {

					//d3.selectAll(".circle")

					// console.log(linegraph_x_domain)
					// console.log(linegraph_y_domain)

					for (var i = 0; i < data.length; i ++) {

						//console.log(data[i].val)
						//console.log(linegraph_y_domain)
						//console.log(data[i].val > linegraph_y_domain[0] && data[i].val < linegraph_y_domain[1])

						if (data[i].rank < linegraph_x_domain[1] && data[i].rank > linegraph_x_domain[0] && data[i].val > linegraph_y_domain[0] && data[i].val < linegraph_y_domain[1]) { 
							//console.log("true")
							addRow("tbody", data[i])
						}

					}

				}

				linegraph_transition();

				d3.select(".linegraph_brush").call(linegraph_brush.clear());

				// add the on click events for the button
				linegraph_clear_button.on('click', function ()
				{
				    // reset everything
					linegraph_x.domain(linegraph_xdomain);
					linegraph_y.domain(linegraph_ydomain);

					linegraph_transition();

					linegraph_clear_button.remove();
				});

				function linegraph_transition() {
					
					linegraph.select(".x.axis")
						.transition()
						.duration(1000)
						.call(linegraph_xAxis);

					linegraph.select(".y.axis")
						.transition()
						.duration(1000)
						.call(linegraph_yAxis);

					linegraph.selectAll("circle")
						.transition()
						.duration(1000)
						.attr("cx", function(d) {
							return linegraph_x(d.rank);
						})
						.attr("cy", function(d) {
							return linegraph_y(d.val);
						});

				}
			}
		}
		
	});

	d3.select("#loading-mask")
	    .style("visibility", "hidden");

}

function deleteRow(tableID) {
	var table = document.getElementById(tableID);
	var rowCount = table.rows.length;

	for (var i = 0; i < rowCount; i++) {
		var row = table.rows[i];
		var checkbox = row.cells[0].childNodes[0];

		if (null!= checkbox && checkbox.checked == true) {
			table.deleteRow(i);
			rowCount--;
			i--;
		}
	}

	$(".tablesorter").trigger("update");
}	

function clearTable(tableID) {
	var table = document.getElementById(tableID);
	var rowCount = table.rows.length;

	for (var i = 0; i < rowCount; i++) {
		table.deleteRow(i);
		rowCount--;
		i--;
	}
}

//export to excel file on click of export button
$(document).ready(function () {

    function exportTableToCSV($table, filename) {

    	// console.log($table)

        var $rows = $table.find('tr:has(td)');
        // console.log($rows)

            // Temporary delimiter characters unlikely to be typed by keyboard
            // This is to avoid accidentally splitting the actual contents
            tmpColDelim = String.fromCharCode(11); // vertical tab character
            tmpRowDelim = String.fromCharCode(0); // null character

            // actual delimiter characters for CSV format
            colDelim = ',';
            rowDelim = '"\r\n"';

            // Grab text from table into CSV formatted string
            csv = '"' + $rows.map(function (i, row) {
                var $row = $(row),
                    $cols = $row.find('td');

                return $cols.map(function (j, col) {
                    var $col = $(col),
                        text = $col.text();

                    return text.replace('"', '""'); // escape double quotes

                }).get().join(tmpColDelim);

            }).get().join(tmpRowDelim)
                .split(tmpRowDelim).join(rowDelim)
                .split(tmpColDelim).join(colDelim) + '"';

            // Data URI
            csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

        $(this)
            .attr({
            'download': filename,
                'href': csvData,
                'target': '_blank'
        });
    }

    // This must be a hyperlink
    $("#export_button").on('click', function (event) {
        // CSV
        exportTableToCSV.apply(this, [$('#table_div>table'), 'export.csv']);

    });	
});

screenshotPreview = function(filename){    

    $("#linegraph").append('<div class="pdf_image"> <object data=' + filename + ' type="application/pdf" width="800px" height="300"> alt : <a href='+ filename + '>test.pdf</a> </object> </div>' )                                                                               
};


$(".dropdown-menu li a").click(function () {
	file = $(this).text()
	output_name = file.split("_HOCKEY");

		d3.csv("/Documents/Bradner_work/hockey-sticks/" + output_name[0] + "_CRC.csv", function(links) {

			var nodes = {};

			var node_list = [];

			var node_names = [];

			// Compute the distinct nodes from the links.
			links.forEach(function(link) {

				link.source = nodes[link.source] || (nodes[link.source] = {name: link.source, fixed: true, weight: +link.weight});
			  	link.target = nodes[link.target] || (nodes[link.target] = {name: link.target, fixed: true, weight: +link.weight});

			  	if (link.node_list) {
			  		node_list.push({name: link.node_list, degree: +link.degree, clique: +link.clique_percentage, enhancerRank: +link.enhancerRank});
			  	}

			});

			var node_length = node_list.length;

			while (node_length--) {
				if (node_list[node_length].degree < 120) {
					node_list.splice(node_length,1)
				}
				else if (node_list[node_length].clique == 0) {
					node_list.splice(node_length,1)
				}
				else {
					node_names.push(node_list[node_length].name)
				}
			}

			var width = 400,
			    margin = 50;

			var radius = 150;

			var svg = d3.select("#crcgraph").append("svg")
			    .attr("width", width)
			    .attr("height", width + 50)
			    .attr("class", "crc_svg")
			    .append("g")
			    .attr("transform", "translate(" + (width/2 + margin)  + "," + (width / 2 + margin + 50) + ")");

			var positions = [];

			var node_opacity_scale = d3.scale.linear()
				.domain(d3.extent(node_list, function(d) { return +d.degree; }))
				.range([1,1])

			var node_size_scale = d3.scale.linear()
				.domain(d3.extent(node_list, function(d) { return +d.enhancerRank; }))
				.range([15, 3])

		    var node_color_scale = d3.scale.linear()
		    	.range(["white", "red"])
		    	.domain([0, 1])


			var node = svg.selectAll("circle")
			    .data(node_list)
			  .enter().append("circle")
			  	.attr("class", "node")
			  	//.style("stroke-width", "20px")
			  	//.style("stroke", "white")
			    .attr("transform", function(d, i) {

			    	positions.push({name: d.name, x: Math.cos(2*Math.PI/node_list.length*i)*radius - 50, y: Math.sin(2*Math.PI/node_list.length*i)*radius - 50, degree: +d.degree});

				  	return "translate(" + (Math.cos(2*Math.PI/node_list.length*i)*radius - 50) +
		 			         "," + (Math.sin(2*Math.PI/node_list.length*i)*radius - 50) + ")";
				})
				.attr("fill", function(d) {
					//console.log(d.clique_percentage)
					return node_color_scale(+d.clique)
				})
				.attr("id", function(d) {
					return d.name;
				})
				.attr("r", function(d) {
					//console.log(d)
					if (d.enhancerRank == 0) {
						return 0;
					}

					else {
						return node_size_scale(d.enhancerRank);
					}
				})
				.attr("opacity", function(d) {

					return node_opacity_scale(d.degree);

				});

			//add whitespace around each circle
			var spacer = svg.selectAll(".spacer")
				.data(node_list)
				.enter().append("circle")
				.attr("class", "spacer")
				.attr("transform", function(d,i) {
					return "translate(" + (Math.cos(2*Math.PI/node_list.length*i)*radius - 50) +
					         "," + (Math.sin(2*Math.PI/node_list.length*i)*radius - 50) + ")";
				})
				.attr("fill", "white")
				.attr("opacity", 1)
				.attr("r", function(d) {
					if (d.degree < 80) {
						return 5;
					}
					else if (d.enhancerRank == 0) {
						return 10;
					}
					else {
						return node_size_scale(d.enhancerRank) + 5
					}
				});

			d3.selection.prototype.moveToFront = function() {
			  return this.each(function(){
			    this.parentNode.appendChild(this);
			  });
			};

			var link_length = links.length;

			while (link_length--) {

				if (node_names.indexOf(links[link_length].target.name) == -1) {
					links.splice(link_length, 1)
					//console.log("hello")
				}

				else if (node_names.indexOf(links[link_length].source.name) == -1) {
					links.splice(link_length, 1)
					//console.log("hello")
				}

				// else {
				// 	console.log("hello")
				// }
			}


			var link_length_2 = links.length;

			//console.log(link_length_2)

			//console.log(links)

			for (var i = 0; i < link_length_2; i++) {
				//console.log(links[i])

				for (var j = 0; j < positions.length; j++) {

					if (links[i].source.name == positions[j].name) {
						links[i].source.x = positions[j].x
						links[i].source.y = positions[j].y
					}
					if (links[i].target.name == positions[j].name) {
						links[i].target.x = positions[j].x
						links[i].target.y = positions[j].y
					}
				}

				//if (i == 1) {console.log(links[i])}
			}

			//console.log(links)

			var stroke_gradient = svg.append("svg:defs")
		    	.append("svg:linearGradient")
		    	.attr("id", "stroke_gradient")
		    	.attr("x1", "0%")
		    	.attr("y1", "0%")
		    	.attr("x2", "100%")
		    	.attr("y2", "0%")
		    	.attr("spreadMethod", "pad");

		    stroke_gradient.append("svg:stop")
		    	.attr("offset", "0%")
		    	.attr("stop-color", "black")
		    	.attr("stop-opacity", 1);

		    stroke_gradient.append("svg:stop")
		    	.attr("offset", "20%")
		    	.attr("stop-color", "black")
		    	.attr("stop-opacity", .2);

		    stroke_gradient.append("svg:stop")
		    	.attr("offset", "80%")
		    	.attr("stop-color", "black")
		    	.attr("stop-opacity", .2);

		    stroke_gradient.append("svg:stop")
		    	.attr("offset", "100%")
		    	.attr("stop-color", "black")
		    	.attr("stop-opacity", .05);

			var stroke_width_scale = d3.scale.linear()
				.range([3, 15])
				.domain(d3.extent(links, function(d) { return +d.weight; }));

			var stroke_opacity_scale = d3.scale.linear()
				.range([0, 1])
				.domain([30, d3.max(links, function(d) { return +d.weight; })]);

			svg.selectAll(".link_line")
				.data(links)
			  .enter()
				.append("path")
			    .attr("class", "link_line")
			    .attr("fill", "url(#stroke_gradient)")
			    .attr("id", function(i, d) { return "link_line" + d; } )
			    .attr("d", function(d){ return drawCurve(d); })
			    .attr("opacity", function(d) {

			    	if (d.weight < 15) {
			    		return 0;
			    	} 
			    	else { 
			    		return 1; 
			    	}

			    });

		    function drawCurve(d) {

		    	var d3LineLinear = d3.svg.line().interpolate("linear");

			    var slope = Math.atan2((d.target.y - d.source.y), (d.target.x - d.source.x));
			    var slopePlus90 = Math.atan2((d.target.y - d.source.y), (d.target.x - d.source.x)) + (Math.PI/2);

			    var sourceX = d.source.x;
			    var sourceY = d.source.y;
			    var targetX = d.target.x;
			    var targetY = d.target.y;

			    var points = [];

			    var bothDirections = -1;
			    $.each(links, function(index, value) {
			      if(value.source == d.target && value.target == d.source) {
			        bothDirections = index;
			        return false;
			      }
			    });

			    if(bothDirections >= 0) {
			      	points.push([sourceX - 0.001*radius*Math.cos(slope), sourceY - 0.001*radius*Math.sin(slope)]);
			      	points.push([targetX + radius*Math.cos(slope) + (stroke_width_scale(links[bothDirections].weight)) * Math.cos(slopePlus90), targetY + radius*Math.sin(slope) + (stroke_width_scale(links[bothDirections].weight)) * Math.sin(slopePlus90)]);
			      	points.push([sourceX - 0.001*radius*Math.cos(slope) + 2 * stroke_width_scale(links[bothDirections].weight) * Math.cos(slopePlus90), sourceY - 0.001*radius*Math.sin(slope) + 2 * stroke_width_scale(links[bothDirections].weight) * Math.sin(slopePlus90)]);
			    } 

			    else {
			      	points.push([sourceX - radius*Math.cos(slope) - stroke_width_scale(d.weight) * Math.cos(slopePlus90), sourceY - radius*Math.sin(slope) - stroke_width_scale(d.weight) * Math.sin(slopePlus90)]);
			      	points.push([targetX  + radius * Math.cos(slope), targetY + radius * Math.sin(slope)]);
			      	points.push([sourceX - radius*Math.cos(slope) + stroke_width_scale(d.weight) * Math.cos(slopePlus90), sourceY - radius*Math.sin(slope) + stroke_width_scale(d.weight) * Math.sin(slopePlus90)]);
			    }

			    return d3LineLinear(points) //+ "Z";
		  	}

			var arc = d3.svg.arc()
				.innerRadius(150)
				.outerRadius(380)
				.startAngle(0)
				.endAngle(360);

			svg.append("path")
				.attr("d", arc)
			    .attr("transform", "translate(-45, -47)")
			    .attr("fill", "white")
			    .attr("class", "circleclear");

			d3.selectAll(".spacer").moveToFront();
			d3.selectAll(".node").moveToFront();

			function computeTextRotation(d) {
			  var angle = d.x - Math.PI / 2;
			  return angle / Math.PI * 180;
			}

			var text = svg.selectAll(".node_text")
				.data(positions)
				.enter().append("text")
		      	//.attr("dy", ".35em")
				.attr("x", function(d) {

					//return d.x - 18

					if (d.x > 0) {
						return d.x + 8;
					} 

					if (d.x < 0) {
						return d.x - 30;
					}


					//return d.x - 18
				})
				.attr("y", function(d) {

					if (d.y > 0) {
						return d.y + 8;
					}

					if (d.y < 0) {
						return d.y - 8;
					}

					//return d.y + 5
				})
				.text(function(d) {
					return d.name;
				})
				.attr("class", "node_text")
				.attr("fill", "black")
				.style("font-size", "12px");


			var gradient = svg.append("svg:defs")
		    	.append("svg:linearGradient")
		    	.attr("id", "gradient")
		    	.attr("x1", "0%")
		    	.attr("y1", "0%")
		    	.attr("x2", "100%")
		    	.attr("y2", "10%")
		    	.attr("spreadMethod", "pad");

		    gradient.append("svg:stop")
		    	.attr("offset", "0%")
		    	.attr("stop-color", "white")
		    	.attr("stop-opacity", 1);

		    gradient.append("svg:stop")
		    	.attr("offset", "100%")
		    	.attr("stop-color", "red")
		    	.attr("stop-opacity", 1);

			svg.append("text")
				.attr("x", width/3 - 170)
				.attr("y", -270)
				.attr("font-size", "16px")
				.attr("font-weight", "bold")
				.attr("text-anchor", "middle")
				.text("Super-enhancer and TF Interaction Network");

		    svg.append("svg:rect")
		    	.attr("width", 150)
		    	.attr("height", 30)
		    	.attr("x", -40)
		    	.attr("y", -260)
		        .attr("class", "legend_gradient")
		    	.style("fill", "url(#gradient)")
		    	.style("stroke-width", "2px")
		    	.style("stroke", "white");

		    svg.append("text")
		    	.attr("x", -25)
		    	.attr("y", -240)
		    	.style("font-weight", "bold")
		    	.text("Clique Percentage")

		    svg.append("text")
		    	.attr("x", -60)
		    	.attr("y", -240)
		    	.text("0%")

		    svg.append("text")
		    	.attr("x", 110)
		    	.attr("y", -240)
		    	.text("100%");

		});
	//d3.select("").remove();

});






