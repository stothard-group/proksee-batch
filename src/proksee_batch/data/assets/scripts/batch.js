////////////////////////////////////////////////////////////////////////////////
// Initial setup
////////////////////////////////////////////////////////////////////////////////

// NOTE: tableData (global variable) is the main data object that contains all the genome data
// Internally we use 'key' instead of 'code_name' to refer to the genome
const genomeData = [...tableData.genomes];
genomeData.forEach( g => g.key = g.code_name );

updateRunInfo();

// Current Sort Settings
let currentSort = { column: null, isDescending: false };
// This will hold a filtered version of the table data when a search is performed
let filteredData = [];
// Current search terms
let searchTerms = [];

// Create viewer
const cgv = new CGV.Viewer('#my-viewer', {height: 500});
window.cgv = cgv;
cgvSetup(cgv);

const initialMessage = 'Click on a genome to view map...';
// Reset map and search
autoResizeMyViewer();
clearMap();
updateSearchCount();

// Initial call to generate table
generateGenomeList();

// Holds currently selected genome key
let selectedKey;

////////////////////////////////////////////////////////////////////////////////
// Main function to generate the genome list
////////////////////////////////////////////////////////////////////////////////

function generateGenomeList() {
    let genomes = (searchTerms.length > 0) ? filteredData : genomeData;
    genomes = sortData(genomes);

    console.log('Generating genome list:', genomes);
    const genomeList = document.getElementById('genome-list-container');
    genomeList.innerHTML = '';
    // const items = Array.isArray(data) ? data : Object.entries(data);
    genomes.forEach((genome, index) => {
        const contigCount = genome.num_contigs;
        const genomeName = highlightSearchTerms(genome.name, searchTerms);
        const genomeDescription = highlightSearchTerms(genome.description, searchTerms);
        const selected = (genome.key == selectedMapKey()) ? 'selected' : '';
        const fileTags = getFileTags(genome);
        genomeList.innerHTML += `
            <div class='genome-item ${selected}' data-key='${genome.key}'>
                <div class='genome-row'>
                    <div class='genome-left'>
                        <div class='genome-name'>${genomeName}</div>
                        <div class='genome-description'>${genomeDescription}</div>
                        <div class='genome-tags'>${fileTags}</div>
                    </div>
                    <div class='genome-right'>
                        <div class='genome-size'>${genome.total_size.toLocaleString()} bp</div>
                        <div class='genome-count'><span class='genome-count-number'>${contigCount}</span> ${(contigCount > 1) ? 'contigs' : 'contig'}</div>
                        <div class='genome-gc'>${(genome.gc_content * 100).toFixed(2)}% GC</div>
                    </div>
                </div>
                <div class='genome-track-listing hidden'>

                </div>
            </div>
        `;
    });

    // Add event listeners to each genome item
    const genomeItems = document.querySelectorAll('.genome-item');
    genomeItems.forEach(item => {
        item.addEventListener('click', function(event) {
            // Could also use event.target.dataset.key
            const key = event.currentTarget.getAttribute('data-key');
            if (key) {
                loadDataForKey(key);
            }
        });
    });
}

function getFileTags(genomeItem) {
    const tagMap = {
        'genbank': 'GenBank',
        'json': 'JSON',
        'blast': 'BLAST',
        'bed': 'BED',
        'vcf': 'VCF',
        'gff': 'GFF',
    };
    console.log("Item:", genomeItem)
    if (!genomeItem.files) return '';
    const fileTags = Object.keys(genomeItem.files).map(fileType => {
        const tagType = tagMap[fileType] || fileType;
        const fileCount = genomeItem.files[fileType].length;
        const filesString = genomeItem.files[fileType].join(', ');
        const tagSuffix = (fileCount > 1) ? `<span class='tag-suffix'> x${fileCount}</span>` : '';
        return (fileCount > 0) ? `<div class='tag tag-${fileType}' title='${filesString}'>${tagType}${tagSuffix}</div>` : '';
    });
    return fileTags.join('');
};

////////////////////////////////////////////////////////////////////////////////
// Sorting functions
////////////////////////////////////////////////////////////////////////////////

const sortButtons = document.querySelectorAll('.sort-btn');
const sortMap = {
    'sort-name': 'name',
    'sort-length': 'total_size',
    'sort-contigs': 'num_contigs',
    'sort-gc': 'gc_content',
};
sortButtons.forEach(btn => {
    btn.addEventListener('click', function(event) {
        const target = event.target;
        const column = sortMap[target.id];
        updateSortSettings(column);
        adjustSortIndicators(target);
        generateGenomeList()
    });
});

function adjustSortIndicators(target) {
    sortButtons.forEach(btn => {
        btn.classList.remove('sort-asc');
        btn.classList.remove('sort-desc');
    });
    if (target.id === 'sort-off') return;
    if (currentSort.isDescending) {
        target.classList.add('sort-desc');
    } else {
        target.classList.add('sort-asc');
    }
}

function updateSortSettings(column) {
    const modifier = currentSort.isDescending ? -1 : 1;

    if (currentSort.column !== column) currentSort.isDescending = false;
    else currentSort.isDescending = !currentSort.isDescending;

    currentSort.column = column;
}

function sortData(data) {
    if (!currentSort.column) return data;
    const column = currentSort.column;
    const modifier = currentSort.isDescending ? -1 : 1;

    const isNumeric = ['total_size', 'num_contigs', 'gc_content'].includes(column);
    let dataArray = [...data];
    dataArray.sort((a, b) => {
        const valueA = isNumeric ? a[column] : a[column].toUpperCase();
        const valueB = isNumeric ? b[column] : b[column].toUpperCase();
        return (valueA < valueB ? -1 : valueA > valueB ? 1 : 0) * modifier;
    });

    const sortedData = [...dataArray];

    return sortedData;
}


////////////////////////////////////////////////////////////////////////////////
// Code for generating Proksee project links
////////////////////////////////////////////////////////////////////////////////

onClick('btn-open-in-proksee', function() {
    const btn = document.getElementById('btn-open-in-proksee');
    if (selectedKey) {
        const genome = genomeData.find(g => g.key == selectedKey);
        if (genome.prokseeUrl) {
            window.open(genome.prokseeUrl, '_blank');
        } else {
            generateProkseeLink(btn, selectedKey)
        }
    }
})

function updateProkseeButton(genome) {
    const btn = document.getElementById('btn-open-in-proksee');
    if (selectedKey) {
        btn.classList.remove('disabled');
        if (genome.prokseeUrl) {
            btn.textContent = 'Open Proksee Project';
            btn.classList.add('link-exists');
        } else {
            btn.textContent = 'Create Proksee Project';
            btn.classList.remove('link-exists');
        }
    } else {
        btn.classList.add('disabled');
    }

}

function generateProkseeLink(element, sampleId) {
    loadScript(`data/genome_maps/${sampleId}.js`, function() {
        if (typeof window.json === 'undefined') {
            console.error('No data found for sample ID:', sampleId);
            alert('Failed to generate Proksee project: No data available for this sample.');
            return;
        }

        const data = { origin: 'proksee-batch', data: JSON.stringify(window.json) };
        const url = 'https://proksee.ca/api/v1/projects.json';
        fetch(url, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data)
        })
        .then(response => response.json())
        .then(data => {
            if (data?.status === 'success' && data?.url) {
                console.log(data)
                const genome = genomeData.find(g => g.key == sampleId);
                genome.prokseeUrl = data.url;
                updateProkseeButton(genome);
            } else {
                console.error(`Failed to create Proksee project: ${data?.error}`);
                alert(`Failed to create Proksee project: ${data?.error}`);
            }
        })
        .catch(error => {
                console.error('Error:', error);
        });
    });
}

function loadScript(scriptUrl, callback) {
    const existingScript = document.querySelector(`script[src="${scriptUrl}"]`);
    if (existingScript) {
        existingScript.remove();
    }
    const script = document.createElement('script');
    script.src = scriptUrl;
    script.type = 'text/javascript';
    document.body.appendChild(script);
    script.onload = function() {
        console.log(`Script loaded: ${scriptUrl}`);
        if (typeof callback === "function") {
            callback();
        }
    };
    script.onerror = function() {
        console.error(`Error loading script: ${scriptUrl}`);
        alert(`Failed to load data for script: ${scriptUrl}`);
    };
}


////////////////////////////////////////////////////////////////////////////////
// Code for loading CGView map from relevant .js files when table rows are clicked.
////////////////////////////////////////////////////////////////////////////////

function loadDataForKey(key) {
    console.log('Clicked on genome:', key);
    selectedKey = key;
    const myViewer = document.querySelector('#my-viewer');
    const mapName = document.querySelector('.map-genome-name');
    showMessage('<div class="spinner-container"><div class="spinner"></div>Loading...</div>');
    const genome = genomeData.find(g => g.key == key);
    mapName.textContent = genome?.name
    highlightSelectedGenome(key);
    hideTrackListing();
    updateProkseeButton(genome);
    var script = document.createElement('script');
    const dataPath = `data/genome_maps/${key}.js`;
    script.src = dataPath;
    script.onload = function() {
        console.log('Loading data file:', dataPath)
        setTimeout(() => {
            myViewer.dataset.key = key;
            cgv.io.loadJSON(window.json);
            // NOTE: resizing fixes an issue where the previous legend is not removed
            addTrackListing(cgv);
            cgv.resize();
            cgv.draw();
            hideMessage();
        }, 1);
    };
    document.head.appendChild(script);
}

function selectedMapKey() {
    const myViewer = document.querySelector('#my-viewer');
    return myViewer.dataset.key;
}

function highlightSelectedGenome(key) {
    document.querySelectorAll('.genome-item').forEach(item => {
        item.classList.remove('selected');
    });
    document.querySelector(`.genome-item[data-key="${key}"]`).classList.add('selected');
}


////////////////////////////////////////////////////////////////////////////////
// Search Bar
////////////////////////////////////////////////////////////////////////////////

const searchInput = document.getElementById('search-input');
searchInput.addEventListener('input', function() {
    const searchCancel = document.getElementById('search-cancel');
    searchCancel.style.display = (searchInput.value !== '') ? 'block' : 'none';
    searchGenomes(searchInput.value);
    updateSearchCount()

});
const searchCancel = document.getElementById('search-cancel');
searchCancel.addEventListener('click', function() {
    searchInput.value = '';
    searchInput.focus();
    searchCancel.style.display = 'none';
    searchGenomes();
    updateSearchCount()
});

function searchGenomes(searchString = '') {
    searchTerms = searchString.trim().toLowerCase().split(/\s+/);
    const dataArray = genomeData;
    console.log('Search terms:', searchTerms);
    if (searchTerms.length > 0) {
        const filteredArray = dataArray.filter(item => {
            itemSearchSpace = item.name + ' ' + item.description;
            return searchTerms.every(term => itemSearchSpace.toLowerCase().includes(term));
        });
        filteredData = filteredArray;
    } else {
        filteredData = genomeData;
    }
    generateGenomeList();
}

function updateSearchCount() {
    const searchInput = document.getElementById('search-input');
    const searchCount = document.getElementById('genome-count');
    totalGenomesCount = genomeData.length;
    filteredGenomesCount = filteredData.length;
    console.log('Total genomes:', totalGenomesCount);
    console.log('Search input:', searchInput.value);
    if (searchInput.value === '') {
        searchCount.textContent = `Genomes total: ${totalGenomesCount}`;
    } else {
        searchCount.textContent = `Genomes found: ${filteredGenomesCount} (of ${totalGenomesCount})`;
    }
}

function highlightSearchTerms(text, terms) {
    if (terms.length === 0) return text;
    const regex = new RegExp(`(${terms.join('|')})`, 'gi');
    return text.replace(regex, '<mark>$1</mark>');
}


////////////////////////////////////////////////////////////////////////////////
// Auto Resize and Drag Bar
////////////////////////////////////////////////////////////////////////////////

function mapResize() {
    const myViewer = document.querySelector('#my-viewer');
    const height = myViewer.clientHeight - 1;
    const width = myViewer.clientWidth - 1;
    // const height = width;
    console.log('Resizing map:', width, height, myViewer.clientHeight);
    cgv.resize(width, height);
}

function autoResizeMyViewer() {
  window.onresize = mapResize;
  window.onload = function () {
    setTimeout( () => {
      mapResize();
      showMessage(initialMessage);
    }, 50);
  }
}

// Drag Bar
const resizeBar = document.querySelector('#resize-bar');
const listPanel = document.querySelector('.section-table');
const mapPanel = document.querySelector('.section-map');
resizeBar.addEventListener('mousedown', function(e) {
    e.preventDefault();
    document.addEventListener('mousemove', onMouseDragMove);
    document.addEventListener('mouseup', onMouseDragUp);
});

function onMouseDragMove(e) {
    const myViewer = document.querySelector('#my-viewer');
    const containerOffsetLeft = listPanel.offsetLeft;
    const newLeftWidth = e.clientX - containerOffsetLeft;
    listPanel.style.width = `${newLeftWidth}px`;
    mapPanel.style.width = `calc(100% - ${newLeftWidth}px)`;
    showMessage('Resizing...');
    myViewer.style.visibility = 'hidden';
}

function onMouseDragUp() {
    const myViewer = document.querySelector('#my-viewer');
    myViewer.style.visibility = 'visible';
    // Only hide message if the map is visible
    const selectedKey = selectedMapKey();
    if (selectedKey) {
        hideMessage();
    } else {
        showMessage(initialMessage);
    }
    mapResize();
    document.removeEventListener('mousemove', onMouseDragMove);
    document.removeEventListener('mouseup', onMouseDragUp);
}

////////////////////////////////////////////////////////////////////////////////
// Run Summary/Details
////////////////////////////////////////////////////////////////////////////////

function updateRunInfo() {
    const inputDirBase = tableData.input_dir.split('/').pop();
    const runName = tableData.run_name || inputDirBase;
    const runDate = tableData.run_date;
    const runVersion = tableData["proksee-batch_version"] || 'Unknown';

    // runDate should be in UTC/GMT format
    // Here we also convert it to local time
    const properDataString = runDate.replace(' ', 'T').replace(' UTC', 'Z');
    // Options to include the local timezone in the output
    let options = {
        timeZoneName: 'short'
    };
    const localRunDate = new Date(properDataString).toLocaleString("en-US", options);

    // Summary
    const runNameEl = document.querySelector('.run-name');
    const runDateEl = document.querySelector('.run-date');
    runNameEl.innerHTML = runName;
    runDateEl.innerHTML = localRunDate;
    // Details
    const detailVersionEl = document.querySelector('.batch-version');
    const detailDateEl = document.querySelector('.batch-run-date');
    const detailInputDirEl = document.querySelector('.batch-input-dir');
    detailVersionEl.innerHTML = runVersion;
    detailDateEl.innerHTML = runDate;
    detailInputDirEl.innerHTML = tableData.input_dir;

    const runSummaryEl = document.querySelector('.run-summary');
    const runDetailsEl = document.querySelector('.run-details');
    runSummaryEl.addEventListener('click', function() {
        if (runSummaryEl.classList.contains('active')) {
            runDetailsEl.classList.remove('active');
            runSummaryEl.classList.remove('active');
        } else {
            runDetailsEl.classList.add('active');
            runSummaryEl.classList.add('active');
        }
    });
}

////////////////////////////////////////////////////////////////////////////////
// Track Listing
////////////////////////////////////////////////////////////////////////////////

function hideTrackListing() {
    document.querySelectorAll('.genome-track-listing').forEach(item => {
        item.classList.add('hidden');
    });
}

function addTrackListing(cgv) {
    const trackListing = document.querySelector(`[data-key='${selectedKey}'] .genome-track-listing`);
    if (trackListing && cgv) {
        trackListing.classList.remove('hidden');
        trackListing.innerHTML = getListing();
    }
}

function getListing() {
  const startOutside = true;
  const separator = "\n";
//   const isLinear = cgv?.format === 'linear';
  const isLinear = false; // Force circular for now
  const label = isLinear ? 'lane' : 'ring'; 
  const collapseTracks = true;
  const direction = isLinear ? (startOutside ? 'top' : 'bottom') : (startOutside ? 'outermost' : 'innermost');
  let text = `<strong>Starting from the ${direction} ${label}:</strong>`;
  text += (separator === '\n') ? "\n" : " ";
  const tracks = cgv?.tracks().filter( (t) => t.visible ) || [];
  // Array of objects with track and slot properties
  // slot is one of: undefined, +, -, -3, -2, -1, +1, +2, +3
  // Tracks (by default) will start from the outside
  const listing = [{track: {name: '<em>Backbone (Contigs)</em>', backbone: true}}];
  for (const track of tracks) {
    // position: 'both', 'inside', 'outside'
    // separateFeaturesBy: 'readingFrame', 'strand', 'none'
    if (track.separateFeaturesBy === 'none' || track.type === 'plot') {
      if (track.position === 'inside') {
        listing.push({track})
      } else {
        listing.unshift({track});
      }
    } else if (track.separateFeaturesBy === 'strand') {
      if (track.position === 'inside') {
        ['+', '-'].forEach( (s) => { listing.push({track, slot: s}) });
      } else if (track.position === 'outside') {
        ['-', '+'].forEach( (s) => { listing.unshift({track, slot: s}) });
      } else {
        listing.unshift({track, slot: '+'});
        listing.push({track, slot: '-'})
      }
    } else if (track.separateFeaturesBy === 'readingFrame') {
      if (track.position === 'inside') {
        ['+3', '+2', '+1', '-1', '-2', '-3'].forEach( (s) => { listing.push({track, slot: s}) });
      } else if (track.position === 'outside') {
        ['-3', '-2', '-1', '+1', '+2', '+3'].forEach( (s) => { listing.unshift({track, slot: s}) });
      } else {
        ['+1', '+2', '+3'].forEach( (s) => { listing.unshift({track, slot: s}) });
        ['-1', '-2', '-3'].forEach( (s) => { listing.push({track, slot: s}) });
      }
    }
  }
  if (!startOutside) {
    listing.reverse();
  }
  let entries = [];
  let slots = [];
  if (collapseTracks) {
    for (let i=0, len=listing.length; i < len; i++) {
      const track = listing[i].track;
      const s = listing[i].slot;
      const next = listing[i+1];
      slots.push(s);
      if (next && next.track === track) {
        continue;
      }
      entries.push(`${displayLabel(i, slots, track.backbone, isLinear)}${track.name}${displaySlots(slots)}`);
      slots = [];
    }
  } else {
    entries = listing.map( (t, i) => `${displayLabel(i, slots, t.track.backbone, isLinear)}${t.track.name}${displaySlots(t.slot)}` );
  }
  return text + entries.join(`${separator}`);
}
// Slots is an array of strings that decribes specific slots
// Returns a empty string or a string cotainin each slot
function displaySlots(slots) {
  let displayText = '';
  if(slots === undefined) return displayText;
  slots = (Array.isArray(slots)) ? slots : [slots];
  if (slots.length > 0 && slots[0] !== undefined) {
    displayText += ` (${slots.join(',')})`;
  }
  return displayText;
}
function displayLabel(index, slots=[], backbone=false, isLinear=false) {
  if (backbone) return '';
  slots = (Array.isArray(slots)) ? slots : [slots];
  let label = isLinear ? 'Lane' : 'Ring';
  if (slots.length > 1) { label += 's'; }
  let numbers = '';
  if (slots.length <= 1) {
    return `${label} ${index+1}: `;
  } else if (slots.length === 2) {
    return `${label} ${index},${index+1}: `;
  } else {
    return `${label} ${index+2 - slots.length}-${index+1}: `;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Helpers: show/hide messages, clear map
////////////////////////////////////////////////////////////////////////////////

function clearMap() {
    cgv.io.loadJSON({cgview: {version: "1.6.0"}}); // Clear the map
    cgv.draw();
    showMessage(initialMessage);
}

function showMessage(message) {
    console.log('SHOWING message');
    const messageContainer = document.querySelector('#my-viewer-message');
    messageContainer.style.opacity = 1;
    messageContainer.style.display = 'flex';
    messageContainer.innerHTML = `<div id='message'>${message}</div>`;
}
function hideMessage() {
    console.log('Hiding message');
    const messageContainer = document.querySelector('#my-viewer-message');
    messageContainer.style.opacity = 0;
    setTimeout(() => {
        messageContainer.style.display = 'none';
    }, 1000); // This should match the transition duration in the CSS (or be a little bigger)
}

////////////////////////////////////////////////////////////////////////////////
// CGView Setup
////////////////////////////////////////////////////////////////////////////////

function cgvSetup(cgv) {
    cgv.on('mousemove.batch', (e) => {
        // const elements = ['caption', 'legendItem'];
        // if (elements.includes(e.elementType)) {
        //   e.element.highlight();
        // }
        if (e.elementType === 'label') {
            const label = e.element;
            label.feature.highlight();
        }
    });
}
