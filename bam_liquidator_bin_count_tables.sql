CREATE DATABASE meta_analysis;
CREATE USER 'counter'@'localhost';

USE meta_analysis;

CREATE TABLE `counts` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `parent_directory` varchar(15) NOT NULL,
  `file_name` varchar(63) NOT NULL,
  `chromosome` varchar(15) NOT NULL,
  `bin` int(11) NOT NULL,
  `count` int(11) NOT NULL,
  `counter_version` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `counter_version` (`counter_version`,`chromosome`,`bin`,`file_name`),
  KEY file (counter_version, file_name)
) ENGINE=MyISAM;

CREATE TABLE errors (id INT NOT NULL AUTO_INCREMENT, file_name VARCHAR(63) NOT NULL,
  chromosome VARCHAR(15) NOT NULL, bin INT NOT NULL, error VARCHAR(512) NOT NULL,
  counter_version INT NOT NULL, create_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP, PRIMARY KEY (id));

CREATE TABLE run (id INT NOT NULL AUTO_INCREMENT, file_name VARCHAR(63) NOT NULL,
                         resume_count INT NOT NULL DEFAULT 0, start_time DATETIME NOT NULL,
                         end_time DATETIME, counter_version INT NOT NULL,
                         PRIMARY KEY (id), UNIQUE KEY (file_name, counter_version));

CREATE TABLE normalized_bins (
  id int NOT NULL AUTO_INCREMENT,
  cell_type varchar(15) NOT NULL,
  chromosome varchar(15) NOT NULL,
  bin int NOT NULL,
  count_fraction double NOT NULL,
  percentile_in_cell_type DOUBLE,
  counter_version int NOT NULL,
  PRIMARY KEY (id),
  UNIQUE KEY plotter (counter_version, chromosome, cell_type, bin),
  KEY overall_view (counter_version, chromosome, bin)
) ENGINE=MyISAM;

GRANT SELECT ON meta_analysis.* TO 'counter'@'localhost';
GRANT SELECT,INSERT ON meta_analysis.counts TO 'counter'@'localhost';
GRANT SELECT,INSERT ON meta_analysis.errors TO 'counter'@'localhost';
GRANT SELECT,INSERT,UPDATE ON meta_analysis.run TO 'counter'@'localhost';
GRANT SELECT,INSERT,UPDATE ON meta_analysis.normalized_bins TO 'counter'@'localhost';
GRANT FILE ON *.* to 'counter'@'localhost';
