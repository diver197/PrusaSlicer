//#undef NDEBUGc
#include <cassert>

#include "PresetBundle.hpp"

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <boost/nowide/cenv.hpp>
#include <boost/nowide/cstdio.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/locale.hpp>

#include <wx/image.h>
#include <wx/choice.h>
#include <wx/bmpcbox.h>
#include <wx/wupdlock.h>

#include "../../libslic3r/libslic3r.h"
#include "../../libslic3r/PlaceholderParser.hpp"
#include "../../libslic3r/Utils.hpp"

namespace Slic3r {

PresetBundle::PresetBundle() :
    prints(Preset::TYPE_PRINT, Preset::print_options()), 
    filaments(Preset::TYPE_FILAMENT, Preset::filament_options()), 
    printers(Preset::TYPE_PRINTER, Preset::printer_options()),
    m_bitmapCompatible(new wxBitmap),
    m_bitmapIncompatible(new wxBitmap)
{
    ::wxInitAllImageHandlers();

    // Create the ID config keys, as they are not part of the Static print config classes.
    this->prints.preset(0).config.opt_string("print_settings_id", true);
    this->filaments.preset(0).config.opt_string("filament_settings_id", true);
    this->printers.preset(0).config.opt_string("print_settings_id", true);
    // Create the "compatible printers" keys, as they are not part of the Static print config classes.
    this->filaments.preset(0).config.optptr("compatible_printers", true);
    this->prints.preset(0).config.optptr("compatible_printers", true);

    this->prints   .load_bitmap_default("cog.png");
    this->filaments.load_bitmap_default("spool.png");
    this->printers .load_bitmap_default("printer_empty.png");

    // FIXME select some icons indicating compatibility.
    this->load_compatible_bitmaps("cog.png", "cog.png");
}

PresetBundle::~PresetBundle()
{
	assert(m_bitmapCompatible != nullptr);
	assert(m_bitmapIncompatible != nullptr);
	delete m_bitmapCompatible;
	m_bitmapCompatible = nullptr;
    delete m_bitmapIncompatible;
	m_bitmapIncompatible = nullptr;
    for (std::pair<const std::string, wxBitmap*> &bitmap : m_mapColorToBitmap)
        delete bitmap.second;
}

void PresetBundle::setup_directories()
{
    boost::filesystem::path dir = boost::filesystem::canonical(Slic3r::data_dir());
    if (! boost::filesystem::is_directory(dir))
        throw std::runtime_error(std::string("datadir does not exist: ") + Slic3r::data_dir());
    std::initializer_list<const char*> names = { "print", "filament", "printer" };
    for (const char *name : names) {
		boost::filesystem::path subdir = (dir / name).make_preferred();
        if (! boost::filesystem::is_directory(subdir) && 
            ! boost::filesystem::create_directory(subdir))
            throw std::runtime_error(std::string("Slic3r was unable to create its data directory at ") + subdir.string());
    }
}

void PresetBundle::load_presets(const std::string &dir_path)
{
    this->prints   .load_presets(dir_path, "print");
    this->filaments.load_presets(dir_path, "filament");
    this->printers .load_presets(dir_path, "printer");
    this->update_multi_material_filament_presets();
}

// Load selections (current print, current filaments, current printer) from config.ini
// This is done just once on application start up.
void PresetBundle::load_selections(const AppConfig &config)
{
    prints.select_preset_by_name(config.get("presets", "print"), true);
    filaments.select_preset_by_name(config.get("presets", "filament"), true);
    printers.select_preset_by_name(config.get("presets", "printer"), true);
    auto   *nozzle_diameter = dynamic_cast<const ConfigOptionFloats*>(printers.get_selected_preset().config.option("nozzle_diameter"));
    size_t  num_extruders   = nozzle_diameter->values.size();   
    this->set_filament_preset(0, filaments.get_selected_preset().name);
    for (unsigned int i = 1; i < (unsigned int)num_extruders; ++ i) {
        char name[64];
        sprintf(name, "filament_%d", i);
        if (! config.has("presets", name))
            break;
        this->set_filament_preset(i, config.get("presets", name));
    }
}

// Export selections (current print, current filaments, current printer) into config.ini
void PresetBundle::export_selections(AppConfig &config)
{
    assert(filament_presets.size() >= 1);
    assert(filament_presets.size() > 1 || filaments.get_selected_preset().name == filament_presets.front());
    config.clear_section("presets");
    config.set("presets", "print",    prints.get_selected_preset().name);
    config.set("presets", "filament", filament_presets.front());
	for (int i = 1; i < filament_presets.size(); ++i) {
        char name[64];
        sprintf(name, "filament_%d", i);
        config.set("presets", name, filament_presets[i]);
    }
    config.set("presets", "printer",  printers.get_selected_preset().name);
}

void PresetBundle::export_selections(PlaceholderParser &pp)
{
    assert(filament_presets.size() >= 1);
    assert(filament_presets.size() > 1 || filaments.get_selected_preset().name == filament_presets.front());
    pp.set("print_preset",    prints.get_selected_preset().name);
    pp.set("filament_preset", filament_presets);
    pp.set("printer_preset",  printers.get_selected_preset().name);
}

bool PresetBundle::load_compatible_bitmaps(const std::string &path_bitmap_compatible, const std::string &path_bitmap_incompatible)
{
    bool loaded_compatible   = m_bitmapCompatible  ->LoadFile(
        wxString::FromUTF8(Slic3r::var(path_bitmap_compatible).c_str()), wxBITMAP_TYPE_PNG);
    bool loaded_incompatible = m_bitmapIncompatible->LoadFile(
        wxString::FromUTF8(Slic3r::var(path_bitmap_incompatible).c_str()), wxBITMAP_TYPE_PNG);
    if (loaded_compatible) {
        prints   .set_bitmap_compatible(m_bitmapCompatible);
        filaments.set_bitmap_compatible(m_bitmapCompatible);
        printers .set_bitmap_compatible(m_bitmapCompatible);
    }
    if (loaded_incompatible) {
        prints   .set_bitmap_compatible(m_bitmapIncompatible);
        filaments.set_bitmap_compatible(m_bitmapIncompatible);
        printers .set_bitmap_compatible(m_bitmapIncompatible);        
    }
    return loaded_compatible && loaded_incompatible;
}

DynamicPrintConfig PresetBundle::full_config() const
{    
    DynamicPrintConfig out;
    out.apply(FullPrintConfig());
    out.apply(this->prints.get_edited_preset().config);
    out.apply(this->printers.get_edited_preset().config);

    auto   *nozzle_diameter = dynamic_cast<const ConfigOptionFloats*>(out.option("nozzle_diameter"));
    size_t  num_extruders   = nozzle_diameter->values.size();

    if (num_extruders <= 1) {
        out.apply(this->filaments.get_edited_preset().config);
    } else {
        // Retrieve filament presets and build a single config object for them.
        // First collect the filament configurations based on the user selection of this->filament_presets.
        std::vector<const DynamicPrintConfig*> filament_configs;
        for (const std::string &filament_preset_name : this->filament_presets)
            filament_configs.emplace_back(&this->filaments.find_preset(filament_preset_name, true)->config);
		while (filament_configs.size() < num_extruders)
            filament_configs.emplace_back(&this->filaments.first_visible().config);
        // Option values to set a ConfigOptionVector from.
        std::vector<const ConfigOption*> filament_opts(num_extruders, nullptr);
        // loop through options and apply them to the resulting config.
        for (const t_config_option_key &key : this->filaments.default_preset().config.keys()) {
            // Get a destination option.
            ConfigOption *opt_dst = out.option(key, false);
            if (opt_dst->is_scalar()) {
                // Get an option, do not create if it does not exist.
                const ConfigOption *opt_src = filament_configs.front()->option(key);
                if (opt_src != nullptr)
                    opt_dst->set(opt_src);
            } else {
                // Setting a vector value from all filament_configs.
                for (size_t i = 0; i < filament_opts.size(); ++ i)
                    filament_opts[i] = filament_configs[i]->option(key);
                static_cast<ConfigOptionVectorBase*>(opt_dst)->set(filament_opts);
            }
        }
    }
    
    static const char *keys[] = { "perimeter", "infill", "solid_infill", "support_material", "support_material_interface" };
    for (size_t i = 0; i < sizeof(keys) / sizeof(keys[0]); ++ i) {
        std::string key = std::string(keys[i]) + "_extruder";
        auto *opt = dynamic_cast<ConfigOptionInt*>(out.option(key, false));
        assert(opt != nullptr);
        opt->value = boost::algorithm::clamp<int>(opt->value, 0, int(num_extruders));
    }

    return out;
}

// Load an external config file containing the print, filament and printer presets.
// Instead of a config file, a G-code may be loaded containing the full set of parameters.
// In the future the configuration will likely be read from an AMF file as well.
// If the file is loaded successfully, its print / filament / printer profiles will be activated.
void PresetBundle::load_config_file(const std::string &path)
{
	if (boost::iends_with(path, ".gcode") || boost::iends_with(path, ".g")) {
		DynamicPrintConfig config;
		config.apply(FullPrintConfig::defaults());
		config.load_from_gcode(path);
		Preset::normalize(config);
		load_config_file_config(path, std::move(config));
		return;
	}

    // 1) Try to load the config file into a boost property tree.
    boost::property_tree::ptree tree;
    try {
        boost::nowide::ifstream ifs(path);
        boost::property_tree::read_ini(ifs, tree);
    } catch (const std::ifstream::failure&) {
        throw std::runtime_error(std::string("The config file cannot be loaded: ") + path);
    } catch (const std::runtime_error&) {
        throw std::runtime_error(std::string("Failed loading the preset file: ") + path);
    }

    // 2) Continue based on the type of the configuration file.
    ConfigFileType config_file_type = guess_config_file_type(tree);
    switch (config_file_type) {
    case CONFIG_FILE_TYPE_UNKNOWN:
        throw std::runtime_error(std::string("Unknown configuration file type: ") + path);   
    case CONFIG_FILE_TYPE_APP_CONFIG:
        throw std::runtime_error(std::string("Invalid configuration file: ") + path + ". This is an application config file.");
	case CONFIG_FILE_TYPE_CONFIG:
	{
		// Initialize a config from full defaults.
		DynamicPrintConfig config;
		config.apply(FullPrintConfig::defaults());
		config.load(tree);
		Preset::normalize(config);
		load_config_file_config(path, std::move(config));
		break;
	}
    case CONFIG_FILE_TYPE_CONFIG_BUNDLE:
		load_config_file_config_bundle(path, tree);
        break;
    }
}

// Load a config file from a boost property_tree. This is a private method called from load_config_file.
void PresetBundle::load_config_file_config(const std::string &path, const DynamicPrintConfig &config)
{
    // 1) Create a name from the file name.
    // Keep the suffix (.ini, .gcode, .amf, .3mf etc) to differentiate it from the normal profiles.
    std::string name = boost::filesystem::path(path).filename().string();

    // 2) If the loading succeeded, split and load the config into print / filament / printer settings.
    // First load the print and printer presets.
    for (size_t i_group = 0; i_group < 2; ++ i_group) {
        PresetCollection &presets = (i_group == 0) ? this->prints : this->printers;
        presets.load_preset(path, name, config).is_external = true;
    }

    // 3) Now load the filaments. If there are multiple filament presets, split them and load them.
    auto   *nozzle_diameter   = dynamic_cast<const ConfigOptionFloats*>(config.option("nozzle_diameter"));
    auto   *filament_diameter = dynamic_cast<const ConfigOptionFloats*>(config.option("filament_diameter"));
    size_t  num_extruders     = std::min(nozzle_diameter->values.size(), filament_diameter->values.size());
    if (num_extruders <= 1) {
        this->filaments.load_preset(path, name, config).is_external = true;
        this->filament_presets.clear();
        this->filament_presets.emplace_back(name);
    } else {
        // Split the filament presets, load each of them separately.
        std::vector<DynamicPrintConfig> configs(num_extruders, this->filaments.default_preset().config);
        // loop through options and scatter them into configs.
        for (const t_config_option_key &key : this->filaments.default_preset().config.keys()) {
            const ConfigOption *other_opt = config.option(key);
            if (other_opt == nullptr)
                continue;
            if (other_opt->is_scalar()) {
                for (size_t i = 0; i < configs.size(); ++ i)
                    configs[i].option(key, false)->set(other_opt);
            } else {
                for (size_t i = 0; i < configs.size(); ++ i)
                    static_cast<ConfigOptionVectorBase*>(configs[i].option(key, false))->set_at(other_opt, 0, i);
            }
        }
        // Load the configs into this->filaments and make them active.
        filament_presets.clear();
        for (size_t i = 0; i < configs.size(); ++ i) {
            char suffix[64];
            if (i == 0)
                suffix[0] = 0;
            else
                sprintf(suffix, " (%d)", i);
            // Load all filament presets, but only select the first one in the preset dialog.
            this->filaments.load_preset(path, name + suffix, std::move(configs[i]), i == 0).is_external = true;
            filament_presets.emplace_back(name + suffix);
        }
    }
}

// Load the active configuration of a config bundle from a boost property_tree. This is a private method called from load_config_file.
void PresetBundle::load_config_file_config_bundle(const std::string &path, const boost::property_tree::ptree &tree)
{
    // 1) Load the config bundle into a temp data.
    PresetBundle tmp_bundle;
    tmp_bundle.load_configbundle(path);

    // 2) Extract active configs from the config bundle, copy them and activate them in this bundle.
    if (tmp_bundle.prints.get_selected_preset().is_default)
        this->prints.select_preset(0);
    else {
        std::string new_name = tmp_bundle.prints.get_selected_preset().name;
        Preset *existing = this->prints.find_preset(new_name, false);
        if (existing == nullptr) {
            // Save under the new_name.
        } else if (existing->config == tmp_bundle.prints.get_selected_preset().config) {
            // Don't save as the config exists in the current bundle and its content is the same.
            new_name.clear();
        } else {
            // Generate a new unique name.
        }
        if (! new_name.empty())
            this->prints.load_preset(path, new_name, std::move(tmp_bundle.prints.get_selected_preset().config));
    }
}

// Load a config bundle file, into presets and store the loaded presets into separate files
// of the local configuration directory.
size_t PresetBundle::load_configbundle(const std::string &path)
{
    // 1) Read the complete config file into a boost::property_tree.
    namespace pt = boost::property_tree;
    pt::ptree tree;
    boost::nowide::ifstream ifs(path);
    pt::read_ini(ifs, tree);

    // 2) Parse the property_tree, extract the active preset names and the profiles, save them into local config files.
    std::vector<std::string> loaded_prints;
    std::vector<std::string> loaded_filaments;
    std::vector<std::string> loaded_printers;
    std::string              active_print;
    std::vector<std::string> active_filaments;
    std::string              active_printer;
    size_t                   presets_loaded = 0;
    for (const auto &section : tree) {
        PresetCollection         *presets = nullptr;
        std::vector<std::string> *loaded  = nullptr;
        std::string               preset_name;
        if (boost::starts_with(section.first, "print:")) {
            presets = &prints;
            loaded  = &loaded_prints;
            preset_name = section.first.substr(6);
        } else if (boost::starts_with(section.first, "filament:")) {
            presets = &filaments;
            loaded  = &loaded_filaments;
            preset_name = section.first.substr(9);
        } else if (boost::starts_with(section.first, "printer:")) {
            presets = &printers;
            loaded  = &loaded_printers;
            preset_name = section.first.substr(8);
        } else if (section.first == "presets") {
            // Load the names of the active presets.
            for (auto &kvp : section.second) {
                if (kvp.first == "print") {
                    active_print = kvp.second.data();
                } else if (boost::starts_with(kvp.first, "filament")) {
                    int idx = 0;
                    if (kvp.first == "filament" || sscanf(kvp.first.c_str(), "filament_%d", &idx) == 1) {
                        if (int(active_filaments.size()) <= idx)
                            active_filaments.resize(idx + 1, std::string());
                        active_filaments[idx] = kvp.second.data();
                    }
                } else if (kvp.first == "printer") {
                    active_printer = kvp.second.data();
                }
            }
        } else if (section.first == "settings") {
            // Load the settings.
            for (auto &kvp : section.second) {
                if (kvp.first == "autocenter") {
                }
            }
        } else
            // Ignore an unknown section.
            continue;
        if (presets != nullptr) {
            // Load the print, filament or printer preset.
            DynamicPrintConfig config(presets->default_preset().config);
            for (auto &kvp : section.second)
                config.set_deserialize(kvp.first, kvp.second.data());
            Preset::normalize(config);
            // Load the preset into the list of presets, save it to disk.
            presets->load_preset(Slic3r::config_path(presets->name(), preset_name), preset_name, std::move(config), false).save();
            ++ presets_loaded;
        }
    }

    // 3) Activate the presets.
    if (! active_print.empty()) 
        prints.select_preset_by_name(active_print, true);
    if (! active_printer.empty())
        printers.select_preset_by_name(active_printer, true);
    // Activate the first filament preset.
    if (! active_filaments.empty() && ! active_filaments.front().empty())
        filaments.select_preset_by_name(active_filaments.front(), true);

    this->update_multi_material_filament_presets();
    for (size_t i = 0; i < std::min(this->filament_presets.size(), active_filaments.size()); ++ i)
        this->filament_presets[i] = filaments.find_preset(active_filaments[i], true)->name;
    return presets_loaded;
}

void PresetBundle::update_multi_material_filament_presets()
{
    // Verify and select the filament presets.
    auto   *nozzle_diameter = static_cast<const ConfigOptionFloats*>(printers.get_edited_preset().config.option("nozzle_diameter"));
    size_t  num_extruders   = nozzle_diameter->values.size();
    // Verify validity of the current filament presets.
    for (size_t i = 0; i < std::min(this->filament_presets.size(), num_extruders); ++ i)
        this->filament_presets[i] = this->filaments.find_preset(this->filament_presets[i], true)->name;
    // Append the rest of filament presets.
//    if (this->filament_presets.size() < num_extruders)
        this->filament_presets.resize(num_extruders, this->filament_presets.empty() ? this->filaments.first_visible().name : this->filament_presets.back());
}

void PresetBundle::export_configbundle(const std::string &path) //, const DynamicPrintConfig &settings
{
    boost::nowide::ofstream c;
    c.open(path, std::ios::out | std::ios::trunc);

    // Put a comment at the first line including the time stamp and Slic3r version.
    c << "# " << Slic3r::header_slic3r_generated() << std::endl;

    // Export the print, filament and printer profiles.
    for (size_t i_group = 0; i_group < 3; ++ i_group) {
        const PresetCollection &presets = (i_group == 0) ? this->prints : (i_group == 1) ? this->filaments : this->printers;
        for (const Preset &preset : presets()) {
            if (preset.is_default || preset.is_external)
                // Only export the common presets, not external files or the default preset.
                continue;
            c << std::endl << "[" << presets.name() << ":" << preset.name << "]" << std::endl;
            for (const std::string &opt_key : preset.config.keys())
                c << opt_key << " = " << preset.config.serialize(opt_key) << std::endl;
        }
    }

    // Export the names of the active presets.
    c << std::endl << "[presets]" << std::endl;
    c << "print = " << this->prints.get_selected_preset().name << std::endl;
    c << "printer = " << this->printers.get_selected_preset().name << std::endl;
    for (size_t i = 0; i < this->filament_presets.size(); ++ i) {
        char suffix[64];
        if (i > 0)
            sprintf(suffix, "_%d", i);
        else
            suffix[0] = 0;
        c << "filament" << suffix << " = " << this->filament_presets[i] << std::endl;
    }

#if 0
    // Export the following setting values from the provided setting repository.
    static const char *settings_keys[] = { "autocenter" };
    c << "[settings]" << std::endl;
    for (size_t i = 0; i < sizeof(settings_keys) / sizeof(settings_keys[0]); ++ i)
        c << settings_keys[i] << " = " << settings.serialize(settings_keys[i]) << std::endl;
#endif

    c.close();
}

// Set the filament preset name. As the name could come from the UI selection box, 
// an optional "(modified)" suffix will be removed from the filament name.
void PresetBundle::set_filament_preset(size_t idx, const std::string &name)
{
    if (idx >= filament_presets.size())
        filament_presets.resize(idx + 1, filaments.default_preset().name);
    filament_presets[idx] = Preset::remove_suffix_modified(name);
}

static inline int hex_digit_to_int(const char c)
{
    return 
        (c >= '0' && c <= '9') ? int(c - '0') : 
        (c >= 'A' && c <= 'F') ? int(c - 'A') + 10 :
        (c >= 'a' && c <= 'f') ? int(c - 'a') + 10 : -1;
}

static inline bool parse_color(const std::string &scolor, unsigned char *rgb_out)
{
    rgb_out[0] = rgb_out[1] = rgb_out[2] = 0;
    const char        *c      = scolor.data() + 1;
    if (scolor.size() != 7 || scolor.front() != '#')
        return false;
    for (size_t i = 0; i < 3; ++ i) {
        int digit1 = hex_digit_to_int(*c ++);
        int digit2 = hex_digit_to_int(*c ++);
        if (digit1 == -1 || digit2 == -1)
            return false;
        rgb_out[i] = (unsigned char)(digit1 * 16 + digit2);
    }
    return true;
}

void PresetBundle::update_platter_filament_ui(unsigned int idx_extruder, wxBitmapComboBox *ui)
{
    if (ui == nullptr)
        return;

    unsigned char rgb[3];
    std::string extruder_color = this->printers.get_edited_preset().config.opt_string("extruder_colour", idx_extruder);
    if (! parse_color(extruder_color, rgb))
        // Extruder color is not defined.
        extruder_color.clear();

    // Fill in the list from scratch.
    ui->Freeze();
    ui->Clear();
    for (size_t i = this->filaments().front().is_visible ? 0 : 1; i < this->filaments().size(); ++ i) {
        const Preset &preset = this->filaments.preset(i);
		if (! preset.is_visible)
			continue;
        bool selected = this->filament_presets[idx_extruder] == preset.name;
		// Assign an extruder color to the selected item if the extruder color is defined.
		std::string   filament_rgb = preset.config.opt_string("filament_colour", 0);
		std::string   extruder_rgb = (selected && !extruder_color.empty()) ? extruder_color : filament_rgb;
		wxBitmap     *bitmap	   = nullptr;
		if (filament_rgb == extruder_rgb) {
			auto it = m_mapColorToBitmap.find(filament_rgb);
			if (it == m_mapColorToBitmap.end()) {
				// Create the bitmap.
				parse_color(filament_rgb, rgb);
				wxImage image(24, 16);
				image.SetRGB(wxRect(0, 0, 24, 16), rgb[0], rgb[1], rgb[2]);
				m_mapColorToBitmap[filament_rgb] = bitmap = new wxBitmap(image);
			} else {
				bitmap = it->second;
			}
		} else {
			std::string bitmap_key = filament_rgb + extruder_rgb;
			auto it = m_mapColorToBitmap.find(bitmap_key);
			if (it == m_mapColorToBitmap.end()) {
				// Create the bitmap.
				wxImage image(24, 16);
				parse_color(extruder_rgb, rgb);
				image.SetRGB(wxRect(0, 0, 16, 16), rgb[0], rgb[1], rgb[2]);
				parse_color(filament_rgb, rgb);
				image.SetRGB(wxRect(16, 0, 8, 16), rgb[0], rgb[1], rgb[2]);
				m_mapColorToBitmap[filament_rgb] = bitmap = new wxBitmap(image);
			} else {
				bitmap = it->second;
			}
		}

		ui->Append(wxString::FromUTF8((preset.name + (preset.is_dirty ? Preset::suffix_modified() : "")).c_str()), (bitmap == 0) ? wxNullBitmap : *bitmap);
        if (selected)
            ui->SetSelection(ui->GetCount() - 1);
    }
    ui->Thaw();
}

// Update the colors preview at the platter extruder combo box.
void PresetBundle::update_platter_filament_ui_colors(unsigned int idx_extruder, wxBitmapComboBox *ui)
{
	this->update_platter_filament_ui(idx_extruder, ui);
	return;

    unsigned char rgb[3];
    std::string extruder_color = this->printers.get_edited_preset().config.opt_string("extruder_colour", idx_extruder);
    if (! parse_color(extruder_color, rgb))
        // Extruder color is not defined.
        extruder_color.clear();

    ui->Freeze();
    for (unsigned int ui_id = 0; ui_id < ui->GetCount(); ++ ui_id) {
        std::string   preset_name        = ui->GetString(ui_id).utf8_str().data();
        size_t        filament_preset_id = size_t(ui->GetClientData(ui_id));
        const Preset *filament_preset    = filaments.find_preset(preset_name, false);
        assert(filament_preset != nullptr);
        // Assign an extruder color to the selected item if the extruder color is defined.
        std::string   filament_rgb       = filament_preset->config.opt_string("filament_colour", 0);
        std::string   extruder_rgb       = (int(ui_id) == ui->GetSelection() && ! extruder_color.empty()) ? extruder_color : filament_rgb;
        wxBitmap     *bitmap             = nullptr;
        if (filament_rgb == extruder_rgb) {
            auto it = m_mapColorToBitmap.find(filament_rgb);
            if (it == m_mapColorToBitmap.end()) {
                // Create the bitmap.
                parse_color(filament_rgb, rgb);
                wxImage image(24, 16);
                image.SetRGB(wxRect(0, 0, 24, 16), rgb[0], rgb[1], rgb[2]);
                m_mapColorToBitmap[filament_rgb] = new wxBitmap(image);
            } else {
                bitmap = it->second;
            }
        } else {
            std::string bitmap_key = filament_rgb + extruder_rgb;
            auto it = m_mapColorToBitmap.find(bitmap_key);
            if (it == m_mapColorToBitmap.end()) {
                // Create the bitmap.
                wxImage image(24, 16);
                parse_color(extruder_rgb, rgb);
                image.SetRGB(wxRect(0, 0, 16, 16), rgb[0], rgb[1], rgb[2]);
                parse_color(filament_rgb, rgb);
                image.SetRGB(wxRect(16, 0, 8, 16), rgb[0], rgb[1], rgb[2]);
                m_mapColorToBitmap[filament_rgb] = new wxBitmap(image);
            } else {
                bitmap = it->second;
            }
        }
        ui->SetItemBitmap(ui_id, *bitmap);
    }
    ui->Thaw();
}

void PresetBundle::set_default_suppressed(bool default_suppressed)
{
    prints.set_default_suppressed(default_suppressed);
    filaments.set_default_suppressed(default_suppressed);
    printers.set_default_suppressed(default_suppressed);
}

} // namespace Slic3r
